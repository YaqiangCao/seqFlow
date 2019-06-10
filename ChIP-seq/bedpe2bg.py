#!/usr/bin/env python
#--coding:utf-8--
"""
bedpe2bg.py
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


def bedpe2model(bg):
    """
    BedGraph format, gzip or not into HTSeq.GenomicArray 
    """
    if bg.endswith(".gz"):
        f = gzip.open(bG, "rb")
    else:
        f = open(bg)
    print(datetime.now(), "Start building model for %s" % bg)
    model = HTSeq.GenomicArray("auto", stranded=False)
    for i, line in enumerate(f):
        if i % 10000 == 0:
            report = "%s lines genome signal read." % i
            cFlush(report)
        line = line.split("\n")[0].split("\t")
        if len(line) < 6:
            continue
        chrom = line[0]
        s = min(int(line[1]), int(line[4]))
        e = max(int(line[2]), int(line[5]))
        iv = HTSeq.GenomicInterval(chrom, s, e)
        model[iv] += 1
    print()
    print(datetime.now(), "Model built for %s" % bg)
    return i, model


def model2bedgraph(t, model, fout):
    with open(fout, "w") as fo:
        for iv, value in model.steps():
            if value > 0:
                value = value / 1.0 / t * 10**6  #RPM
                line = [iv.chrom, iv.start, iv.end, value]
                line = list(map(str, line))
                fo.write("\t".join(line) + "\n")


def bedpe2bdg(f):
    fo = f.split("/")[-1].replace(".bedpe", ".bdg")
    if os.path.isfile(fo):
        return
    t, model = bedpe2model(f)
    model2bedgraph(t, model, fo)


def main():
    fs = glob("../7.seperated/*.bedpe")
    Parallel(n_jobs=len(fs))(delayed(bedpe2bdg)(f) for f in fs)


if __name__ == '__main__':
    start_time = datetime.now()
    main()
    elapsed = datetime.now() - start_time
    fn = os.path.basename(__file__)
    print datetime.now(), "The process is done for %s,time used:%s" % (fn,
                                                                       elapsed)
