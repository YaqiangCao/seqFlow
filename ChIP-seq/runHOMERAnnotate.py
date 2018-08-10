#!/usr/bin/env python2.7
#--coding:utf-8 --
"""
run_HOMER_annotePeak.py
Annotating peaks with HOMER and get peaks-genes relations
2015-10-08: modified.
"""

#
import glob, os, time, subprocess, copy
from datetime import datetime

#
from joblib import Parallel, delayed
import pandas as pd
import numpy as np

#
from Biolibs.rel.General.logger import getlogger
from Biolibs.rel.General.countLines import count_lines
date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
logger = getlogger(fn=os.getcwd() + "/" + date.strip() + "_" +
                   os.path.basename(__file__) + ".log")

__author__ = "CAO Yaqiang"
__date__ = ""
__modified__ = "2015-09-21"
__license__ = "GPL"
__version__ = "0.1"
__email__ = "caoyaqiang0410@gmail.com"


def runHOMER(bed, genome="hg38"):
    name = bed.split("/")[-1].split(".")[0]
    stat = name + ".stat"
    anobed = name + "_ano.bed"
    if os.path.exists(anobed):
        print anobed, "generated, return."
        return
    cmd = "annotatePeaks.pl {bed} {genome} -go {go_dir} > {anobed}".format(
        bed=bed,
        genome=genome,
        stat=stat,
        genome_dir="genome_" + name,
        go_dir="go_" + name,
        anobed=anobed)
    try:
        logger.info(cmd)
        subprocess.call(cmd, shell=True)
    except:
        return


def parseTarget():
    fs = glob.glob("*ano.bed")
    for f in fs:
        print f
        pre = f.split(".")[0]
        mat = pd.read_table(f, index_col=0)
        gs = {}
        for i in mat.index:
            g = mat.loc[i, "Nearest PromoterID"]
            if g not in gs:
                gs[g] = set()
            gs[g].add(i)
        with open(pre + ".gmt", "w") as nf:
            for g, rs in gs.items():
                line = g + "\t" + ",".join(rs) + "\n"
                nf.write(line)


def main():
    beds = glob.glob("../../3.Modules/*.bed")
    map(runHOMER, beds)
    parseTarget()


if __name__ == "__main__":
    startTime = datetime.now()
    main()
    elapsed = datetime.now() - startTime
    print "The process is done."
    print "Time used:", elapsed
