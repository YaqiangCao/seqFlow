#!/usr/bin/env python2.7
#--coding:utf-8 --
"""
"""

__author__ = "CAO Yaqiang"
__date__ = "2015-05-18"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#general library
import glob, os, time, sys, gc

#3rd library
import matplotlib as mpl
mpl.use("pdf")
mpl.rcParams["pdf.fonttype"] = 42
import HTSeq
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import pylab
import brewer2mpl
import seaborn as sns
colors = brewer2mpl.get_map('Set3', 'qualitative', 8).mpl_colors

#my own library
from Biolibs.rel.General.logger import getlogger

#global settings
#logger
date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
logger = getlogger(fn=os.getcwd() + "/" + date.strip() + "_" +
                   os.path.basename(__file__) + ".log")


def flush(i):
    if i % 1000 == 0:
        report = "\r%dk reads parsed" % (i / 1000)
        sys.stdout.write(report)
        sys.stdout.flush()


def adjPets(d):
    #chr_1 == chr_2
    if d[0] == d[3]:
        if d[1] + d[2] > d[4] + d[5]:
            d = [d[0], d[4], d[5], d[3], d[1], d[2], d[6], d[7], d[9], d[8]]
            return d
    #chr_1 > chr_2
    if d[0] > d[3]:
        d = [d[3], d[4], d[5], d[0], d[1], d[2], d[6], d[7], d[9], d[8]]
        return d
    return d


def pairReads(bed):
    sample = bed.split(".")[0]
    fn = sample + ".bedpe"
    if os.path.exists(fn):
        print "%s exists" % fn
        return
    logger.info("PAIRING MAPPED READS. Reading %s" % bed)
    data = {}
    for i, line in enumerate(open(bed)):
        flush(i)
        line = line.split("\n")[0].split("\t")
        rid = int(line[3].split("/")[0].split(".")[1])
        #columns = ["chr_1","start_1","end_1","chr_2","start_2","end_2","rid","-", "strand_1","strand_2" ]
        if rid not in data:
            nr = [
                line[0],
                int(line[1]),
                int(line[2]), False, False, False, rid, "-", line[5], False
            ]
            data[rid] = nr
        else:
            data[rid][3] = line[0]
            data[rid][4] = int(line[1])
            data[rid][5] = int(line[2])
            data[rid][9] = line[5]
    with open(fn, "w") as f:
        for i, rid in enumerate(data.keys()):
            if False in data[rid]:
                continue
            r = adjPets(data[rid])
            line = "\t".join(map(str, r)) + "\n"
            f.write(line)
            if i % 1000 == 0:
                report = "\r%dk pairs writen" % (i / 1000)
                sys.stdout.write(report)
                sys.stdout.flush()
    del data
    print


def main():
    bed = "K562_GSM1551623.bed"
    pairReads(bed)


if __name__ == "__main__":
    main()
