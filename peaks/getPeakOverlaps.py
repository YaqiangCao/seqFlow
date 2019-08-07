#!/usr/bin/env python3.6
#--coding:utf-8 --

"""
getPeakOverlaps.py
2019-08-07:
"""

import os,subprocess
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


def getJI(fi,fj):
    """
    Get the overlap stats
    """
    ci = getN(fi)
    cj = getN(fj)
    cmd = "intersectBed -a %s -b %s -u > ./tmp.bed" % ( fi, fj)
    subprocess.call(cmd, shell=True)
    cij = getN("./tmp.bed")
    ji = cij / 1.0 / (ci + cj - cij)
    return ci,cj,cij,ji
    

def evulateReps(fin, fout):
    mat = pd.read_csv(fin, index_col=0, sep="\t")
    ds = {}
    ns = list(mat.index)
    jis = []
    for i in range(len(ns)):
        ni = "_".join(ns[i].split("_")[:-1])
        #if ni == "WT_6_Th2_Neg":
        #    continue
        for j in range(i + 1, len(ns)):
            nj = "_".join(ns[j].split("_")[:-1])
            if ni == nj:
                ds[ni] = mat.loc[ns[i], ns[j]]
            else:
                jis.append( mat.loc[ns[i],ns[j]] )
    ds = pd.Series(ds)
    ss = list(ds.index)
    ss.sort()
    ds = ds[ss]
    fig, ax = pylab.subplots()
    x = range(len(ds))
    y = ds.values
    ax.vlines(x=x, ymin=0, ymax=y, color="skyblue")
    ax.plot(x, y, "o")
    ax.set_ylabel("Jaccard index of overlap peaks")
    ax.set_xticks(x)
    ax.set_xticklabels(ds.index, rotation=45, fontsize=6, ha="right")
    fig.tight_layout()
    ax.set_title("Mean JI for reps: %f \n Mean JI for non-reps: %f " % (ds.mean(),np.mean(jis)))
    pylab.savefig(fout + ".pdf")


def pre(key="H1_H3K4me3"):
    fas = glob("../1.peaks/*%s*.bed"%key)
    fb = glob("../0.ENCODE/*%s*.bed"%key)[0]
    ds = {}
    for f in fas:
        pre = f.split("/")[-1].split(".")[0]
        d = "_".join( pre.split("_")[:2] )
        tool = "_".join( pre.split("_")[2:] )
        ds[tool] = {"f":f,"ref":fb}
    return ds

def getStat( ds,pre ):
    data = {}
    for tool in ds.keys():
        ci,cj,cij,ji = getJI( ds[tool]["f"],ds[tool]["ref"] )
        data[tool] = {
            "Peaks": ci,
            "ENCODE_Peaks": cj,
            "OverlapPeaks": cij,
            "OverlapRatio": cij /1.0/ ci,
            "JaccardIndex": ji,
            }
    data = pd.DataFrame(data).T
    data.to_csv( pre +".txt",sep="\t")

def main():
    key = "H1_H3K4me3"
    ds = pre(key)
    getStat( ds,key )
    key = "K562_H3K4me3"
    ds = pre(key)
    getStat( ds,key )


if __name__ == "__main__":
    main()
