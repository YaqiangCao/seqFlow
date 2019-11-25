#!/usr/bin/env python2.7
#--coding:utf-8--
"""

"""

__author__ = "CAO Yaqiang"
__date__ = ""
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import os, time, gzip, sys
from glob import glob
from datetime import datetime
from copy import deepcopy

#computating setting
import pandas as pd
import numpy as np
from tqdm import tqdm

#cLoops2 settings
from cLoops2.settings import *


def checkOverlap(ra, rb):
    """
    check the overlap of a region for the same chromosome
    """
    if rb[0] <= ra[0] <= rb[1] or rb[0] <= ra[1] <= rb[1]:
        return True
    if ra[0] <= rb[0] <= ra[1] or ra[0] <= rb[1] <= ra[1]:
        return True
    return False


def getSet(f):
    rs = {}
    for line in open(f):
        line = line.split("\n")[0].split("\t")
        c = line[0]
        s = int(line[1])
        e = int(line[2])
        if c not in rs:
            rs[c] = []
        rs[c].append([s, e])
    #sort the data and remove the small overlapped ones
    for c, ts in rs.items():
        nts = sorted(ts, key=lambda x: x[1] - x[0], reverse=True)
        #nts = sorted(ts,key=lambda x: x[1]-x[0],reverse=False)
        toremove = set()
        ts = []
        for i in range(0, len(nts)):
            if i in toremove:
                continue
            flag = False
            for j in range(i + 1, len(nts)):
                if checkOverlap(nts[i], nts[j]):
                    flag = True
                    toremove.add(j)
            ts.append(nts[i])
        rs[c] = ts
    return rs


def getOverlap(rs, r):
    for i in range(len(rs)):
        if checkOverlap(rs[i], r):
            return i
    return None


def findTADOverlap(rs, f):
    c, o = 0, 0
    for i, line in enumerate(open(f)):
        if i == 0:
            continue
        line = line.split("\n")[0].split("\t")
        c += 1
        chrom = line[1]
        if chrom not in rs:
            continue
        a = [int(line[2]), int(line[3])]
        b = [int(line[5]), int(line[6])]
        pa = getOverlap(rs[chrom], a)
        pb = getOverlap(rs[chrom], b)
        if pa is None or pb is None:
            continue
        if pa == pb:
            o += 1
    return f.split("/")[-1].split("_loops")[0], c, o


def get():
    f = "GM12878_TAD_hg38.bed"
    rs = getSet(f)
    loopfs = glob("../1.loops/*.txt")
    ds = {}
    for f in loopfs:
        n, c, o = findTADOverlap(rs, f)
        ds[n] = {"total": c, "withInTAD": o}
    ds = pd.DataFrame(ds).T
    c = ds["withInTAD"] / ds["total"]
    ds.to_csv("stat.txt", sep="\t")


def plot():
    mat = pd.read_csv("stat.txt", index_col=0, sep="\t")
    c = mat["withInTAD"] / mat["total"]
    c = c.sort_values(ascending=True, inplace=False)
    fig, ax = pylab.subplots(figsize=(3 * 0.8, 2.75 * 0.8))
    sns.barplot(x=list(range(len(c))), y=c, color=colors[1])
    ax.set_ylabel("loops within TAD ratio")
    ax.set_xticklabels(list(c.index), rotation="45", ha="right")
    pylab.savefig("withInTADPETsRatio.pdf")


def main():
    #get()
    plot()


if __name__ == "__main__":
    main()
