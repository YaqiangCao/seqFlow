#!/usr/bin/env python3.6
#--coding:utf-8 --
"""
getPeakOverlaps.py
2019-08-07:
"""

import os, subprocess
from glob import glob

#3rd
import matplotlib as mpl
mpl.use("pdf")
import pylab
import brewer2mpl
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm

#global settings
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["figure.figsize"] = (4, 2.75)
mpl.rcParams["font.size"] = 10.0
sns.set_style("white")
colors = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors


def getN(f):
    c = open(f).read().split("\n")
    c = [t for t in c if t.strip() != ""]
    c = len(c)
    return c


def getJI(fi, fj):
    """
    Get the overlap stats
    """
    ci = getN(fi)
    cj = getN(fj)
    cmd = "intersectBed -a %s -b %s -u > ./tmp.bed" % (fi, fj)
    subprocess.call(cmd, shell=True)
    cij = getN("./tmp.bed")
    ji = cij / 1.0 / (ci + cj - cij)
    return ci, cj, cij, ji



def main():
    data = {}
    for cut in np.arange(0.05,0.5,0.01):
        fa = "./beds/trac1_%.2f_border.bed"%cut
        fb = "./beds/trac2_%.2f_border.bed"%cut
        ci,cj,cij,ji = getJI( fa,fb)
        data[cut] = {
            "trac1":ci,
            "trac2":cj,
            "Overlaps": cij,
            "OverlapRatio": cij / 1.0 / ci,
            "JaccardIndex": ji,
        }
    data = pd.DataFrame(data).T
    data.to_csv("bordersOverlap.txt", sep="\t")



if __name__ == "__main__":
    main()
