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
from cLoops2.settings import *


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
    n = f.split("/")[-2]
    ds = []
    stat = {}
    for line in gzip.open(f,"rt"):
        line = line.split("\n")[0].split("\t")
        s = int(line[1])
        e = int(line[2])
        cid = line[3]
        m = (s + e) / 2
        d = e - s
        ds.append(d)
        if cid not in stat:
            stat[cid] = {
                "total": 0,
                "sP": 0,
                "cN": 0,
                "other": 0,
                "ds": [],
            }
        stat[cid]["total"] += 1
        stat[cid]["ds"].append(d)
        if d <= dfilter[0]:  #subnucleosome-sized particles
            stat[cid]["sP"] += 1
        elif dfilter[1] <= d <= dfilter[2]:  #canonical nucleosomes
            stat[cid]["cN"] += 1
        else:
            stat[cid]["other"] += 1
    for cid in stat.keys():
        stat[cid]["fragmentLengthMean"] = np.mean(stat[cid]["ds"])
        stat[cid]["fragmentLengthStd"] = np.std(stat[cid]["ds"])
        del stat[cid]["ds"]
    plotDis(ds, n)
    stat = pd.DataFrame(stat).T
    stat.to_csv(n + "_stat.txt", sep="\t")


def getStat2(f, dfilter=[137,157]):
    """
    """
    n = f.split("/")[-2]
    ds = []
    stat = {}
    for line in gzip.open(f,"rt"):
        line = line.split("\n")[0].split("\t")
        s = int(line[1])
        e = int(line[2])
        cid = line[3]
        m = (s + e) / 2
        d = e - s
        ds.append(d)
        if cid not in stat:
            stat[cid] = {
                "total": 0,
                "nucleosomes": 0,
                "ds": [],
            }
        stat[cid]["total"] += 1
        stat[cid]["ds"].append(d)
        if dfilter[0] <= d <= dfilter[1]:  #canonical nucleosomes
            stat[cid]["nucleosomes"] += 1
    for cid in stat.keys():
        stat[cid]["fragmentLengthMean"] = np.mean(stat[cid]["ds"])
        stat[cid]["fragmentLengthStd"] = np.std(stat[cid]["ds"])
        stat[cid]["nucleosomeRatio"] = stat[cid]["nucleosomes"] / stat[cid]["total"]
        del stat[cid]["ds"]
    plotDis(ds, n)
    stat = pd.DataFrame(stat).T
    stat.to_csv(n + "_stat2.txt", sep="\t")


def selCells():
    fs = glob("*_stat2.txt")
    for f in fs:
        n = f.split("_stat2")[0]
        mat = pd.read_csv(f,index_col=0,sep="\t")
        sa = mat["total"]
        sa = np.log10(sa)
        sb = mat["nucleosomeRatio"]
        data = pd.DataFrame({"total reads, log10": sa, "137-157 bp ratio": sb})
        fig, ax = pylab.subplots(figsize=(3.2,2.2))
        """
        g = sns.jointplot(data=data,
                          x="total reads, log10",
                          y="137-157 bp ratio",
                          kind="kde",
                          fill=True,
                          marginal_ticks=True)
        """
        #sns.scatterplot(data=data,
        #                  x="total reads, log10",
        #                  y="137-157 bp ratio",ax=ax,size=1)
        ax.scatter(sa,sb,s=1)
        ax.set_xlabel("total reads, log10")
        ax.set_ylabel("137-157 bp ratio")
        ax.set_title(n)
        ax.set_xlim([2,5.5])
        ax.set_ylim([0.12,0.3])
        pylab.savefig(f"{n}_rawCellStat.pdf")

def statWell():
    fs = glob("*_stat2.txt")
    for f in fs:
        n = f.split("_stat2")[0]
        mat = pd.read_csv(f,index_col=0,sep="\t")
        cells = list(mat.index)
        wells = {}
        for cid in cells:
            wid = "_".join(cid.split("_")[:-1])
            if wid not in wells:
                wells[wid] = []
            wells[wid].append( cid )
        cids = []
        for wid, cs in wells.items():
            s = {}
            for cid in cs:
                s[cid] = mat.loc[cid,"nucleosomes"]
            s = pd.Series(s)
            s = s.sort_values(inplace=False, ascending=False)
            cids.extend(list(s.index[:5]))
        print(len(cids))
        with open(n+"_selCells.list","w") as fo:
            fo.write("\n".join(cids))

 
if __name__ == '__main__':
    fs = glob("../2.pre/*/*tsv.gz")
    Parallel(n_jobs=min(len(fs), 10))(delayed(getStat2)(f) for f in fs)
    selCells()
    statWell()
