#!/usr/bin/env python
#--coding:utf-8--
"""
bedpe2bg.py
2019-06-24: click added to use command line interface
"""

__author__ = "CAO Yaqiang"
__date__ = ""
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#sys
import os, gzip, sys, random,time
from datetime import datetime
from glob import glob
from collections import Counter

#3rd library
import click
import HTSeq
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

# this
from utils import cFlush,getLogger,callSys,PET

#global setting
#logger
date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
logger = getLogger(fn=os.getcwd() + "/" + date.strip() + "_" +
                   os.path.basename(__file__) + ".log")


def bedpe2model(bg, mapq=1):
    """
    Convet BEDPE format file into HTSeq.GenomicArray to get the genomic coverage.
    Only non-redundant reads will be kept.

    Parameteres
    ----
    bg: str, .bedpe or .bedpe.gz file
    mapq: int, Bowtie2 MAPQ cutoff to filter reads.

    Returns
    ----
    HTSeq.GenomicArray
    """
    rs = set()
    if bg.endswith(".gz"):
        fh = gzip.open(bg, "rb")
    else:
        fh = open(bg)
    logger.info("Start building model for %s, with MAPQ cutoff >=%s" % (bg,mapq))
    model = HTSeq.GenomicArray("auto", stranded=False)
    t = 0
    for i, line in enumerate(fh):
        if i % 10000 == 0:
            report = "%s lines genome signal read." % i
            cFlush(report)
        line = line.split("\n")[0].split("\t")
        try:
            pet = PET(line) 
        except:
            logger.error("%s from %s is not a BEDPE record"%(line,bg))
        if not pet.cis or "_" in pet.chromA:
            continue
        if pet.mapq < mapq:
            continue
        t += 1
        r = (pet.chromA, pet.mid, pet.mid + 1)
        if r not in rs:
            iv = HTSeq.GenomicInterval(pet.chromA, pet.start, pet.end)
            model[iv] += 1
            rs.add(r)
    logger.info("%s:totalReads:%s;nonRedudant:%s" % (bg, t, len(rs)))
    return len(rs), model


def model2bedgraph(t, model, fout):
    """
    Converting HTSeq.GenomicArray to RPKM measured BEDGRAPH file.

    Parameters
    ---
    t: int, number of total reads for the sample
    model: HTSeq.GenomicArray
    fout: str, output file name
    """
    with open(fout, "w") as fo:
        for iv, value in model.steps():
            if value > 0:
                value = value / 1.0 / t * 10**6  #RPM
                line = [iv.chrom, iv.start, iv.end, value]
                line = list(map(str, line))
                fo.write("\t".join(line) + "\n")


def bedpe2bdg(f,mapq=1):
    """
    Convert BEDPE file to BEDGRAPH file.

    Parameter
    ---
    f: str, input file, .bedpe.gz  
    mapq: int, mapq cutoff
    """
    fo = f.split("/")[-1].replace(".bedpe.gz", ".bdg")
    if os.path.isfile(fo):
        return
    t, model = bedpe2model(f,mapq)
    model2bedgraph(t, model, fo)


@click.command()
@click.option("-dir",required=True,help="Directory for the .bedpe.gz files.")
@click.option("-mapq",default=1,help="MAPQ cutoff for filtering PETs.")
@click.option("-cpu",default=10,help="Number of CPUs to finish the job, default is set to 10.")
def main(dir,cpu,mapq):
    """
    Converting .bedpe.gz files from other directory into this directory .bedgraph files.
    """
    if not os.path.exists(dir):
        logger.error("%s not exists!"%dir)
    fs = glob("%s/*.bedpe.gz"%dir)
    cpu = min(cpu,len(fs))
    Parallel(n_jobs=cpu)(delayed(bedpe2bdg)(f,mapq) for f in fs)


if __name__ == '__main__':
    start_time = datetime.now()
    main()
    elapsed = datetime.now() - start_time
    fn = os.path.basename(__file__)
    print(datetime.now(), "The process is done for %s,time used:%s" % (fn, elapsed))
