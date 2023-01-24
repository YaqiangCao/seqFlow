#!/usr/bin/env python2.7
#--coding:utf-8--
"""
tsvStat.py
scMNase-seq data aggreagtion stat.
"""

__author__ = "CAO Yaqiang"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import os
import gzip
from glob import glob

#3rd library
import numpy as np
import pandas as pd
import seaborn as sns
import pylab
from joblib import Parallel, delayed


def plotDis(ds, name):
    fig, ax = pylab.subplots(figsize=(3.2,2.2))
    sns.kdeplot(ds, ax=ax)
    ax.axvline(x=80, linewidth=1, linestyle="--", color="gray")
    ax.text(80, 0.001, "80bp")
    ax.axvline(x=120, linewidth=1, linestyle="--", color="gray")
    ax.text(120, 0.001, "120bp")
    ax.axvline(x=140, linewidth=1, linestyle="--", color="gray")
    ax.text(140, 0.001, "140bp")
    ax.axvline(x=180, linewidth=1, linestyle="--", color="gray")
    ax.text(180, 0.001, "180bp")
    ax.set_xlim([0, 300])
    ax.set_title("fragment length")
    ax.set_xlabel("length")
    ax.set_ylabel("density")
    pylab.savefig(f"{name}_fragmentLength.pdf")


def getStat(f, dfilter=[120, 140, 180]):
    """
    Get the name, redundancy, total PETs, PETs distance mean, distance std for a bedpe file.
    """
    n = f.split("/")[-1].split(".bedpe")[0]
    ds = []
    stat = {
                "total": 0,
                "sP": 0,
                "cN": 0,
                "other": 0,
                "ds": [],
            }

    for line in gzip.open(f,"rt"):
        line = line.split("\n")[0].split("\t")
        s = int(line[1])
        #e = int(line[2])
        e = int(line[5])
        m = (s + e) / 2
        d = e - s
        ds.append(d)
        stat["total"] += 1
        if d <= dfilter[0]:  #subnucleosome-sized particles
            stat["sP"] += 1
        elif dfilter[1] <= d <= dfilter[2]:  #canonical nucleosomes
            stat["cN"] += 1
        else:
            stat["other"] += 1
    plotDis(ds, n)
    #stat = pd.DataFrame(stat).T
    stat = pd.Series(stat)
    stat.to_csv(n + "_stat.txt", sep="\t")


if __name__ == '__main__':
    fs = glob("*bedpe.gz")
    Parallel(n_jobs=min(len(fs), 10))(delayed(getStat)(f) for f in fs)
