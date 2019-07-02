#!/usr/bin/env python
#--coding:utf-8--
"""
bedpe2bg.py
2019-07-02: extension to 150 bp default added.
"""

__author__ = "CAO Yaqiang"
__date__ = ""
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#sys
import os, gzip, sys, random
from datetime import datetime
from glob import glob
from collections import Counter

#3rd library
import HTSeq
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

# this
from utils import *

#global setting
#logger
date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
logger = getLogger(fn=os.getcwd() + "/" + date.strip() + "_" +
                   os.path.basename(__file__) + ".log")


def bed2model(bg,mapq=1,noRedu=True,ext=150):
    """
    Convet BED format file into HTSeq.GenomicArray to get the genomic coverage.
    Only non-redundant reads will be kept.

    Parameteres
    ----
    bg: str, .bed or .bed.gz file
    mapq: int, Bowtie2 MAPQ cutoff to filter reads.
    noRedu: bool, whether to keep redundant reads

    Returns
    ----
    HTSeq.GenomicArray
 
    BedGraph format, gzip or not into HTSeq.GenomicArray 
    """
    rs = set()
    if bg.endswith(".gz"):
        f = gzip.open(bg, "rb")
    else:
        f = open(bg)
    logger.info("Start building model for %s, with MAPQ cutoff >=%s" %
                (bg, mapq))
    model = HTSeq.GenomicArray("auto", stranded=False)
    t = 0
    for i, line in enumerate(f):
        if i % 10000 == 0:
            report = "%s lines genome signal read." % i
            cFlush(report)
        line = line.split("\n")[0].split("\t")
        if len(line) < 3:
            continue
        try:
            chrom = line[0]
            s = int(line[1])
            e = int(line[2])
        except:
            continue
        if int(line[4]) < mapq:
            continue
        t += 1
        r = (chrom, s, e)
        if noRedu:
            if r not in rs:
                if line[5] == "+":
                    e = s+ext
                else:
                    s = max(0,e - ext)
                iv = HTSeq.GenomicInterval(chrom, s, e)
                model[iv] += 1
                rs.add(r)
        else:
            iv = HTSeq.GenomicInterval(chrom, s, e)
            model[iv] += 1
    print("%s:totalReads:%s;nonRedudant:%s" % (f, i, len(rs)))
    logger.info("%s:totalReads:%s;nonRedudant:%s" % (f, i, len(rs)))
    return len(rs), model


def model2bedgraph(t, model, fout):
    with open(fout, "w") as fo:
        for iv, value in model.steps():
            if value > 0:
                value = value / 1.0 / t * 10**6  #RPM
                line = [iv.chrom, iv.start, iv.end, value]
                line = list(map(str, line))
                fo.write("\t".join(line) + "\n")


def bed2bdg(f):
    fo = f.split("/")[-1].replace(".bed.gz", ".bdg")
    if os.path.isfile(fo):
        return
    t, model = bed2model(f)
    model2bedgraph(t, model, fo)


def main():
    fs = glob("../3.beds/*.bed.gz")
    Parallel(n_jobs=len(fs))(delayed(bed2bdg)(f) for f in fs)


if __name__ == '__main__':
    start_time = datetime.now()
    main()
    elapsed = datetime.now() - start_time
    fn = os.path.basename(__file__)
    print(datetime.now(),
          "The process is done for %s,time used:%s" % (fn, elapsed))
