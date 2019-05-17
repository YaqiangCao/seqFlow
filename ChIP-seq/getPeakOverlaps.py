#!/usr/bin/env python3.6
#--coding:utf-8 --
"""
getPeakOverlaps.py
2019-05-13:
"""

import os, subprocess
from glob import glob

#3rd
import pylab
import brewer2mpl
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm

#global settings
import matplotlib as mpl
mpl.use("pdf")
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


def getJIs(fs, fnOut):
    """
    Get the Jaccard Index between peaks.
    """
    ds = {}
    for fi in tqdm(fs):
        #ni = fi.split("/")[-1].split(".")[0]
        ni = fi.split("/")[-2]
        ds[ni] = {}
        ci = getN(fi)
        for fj in fs:
            #nj = fj.split("/")[-1].split(".")[0]
            nj = fj.split("/")[-2]
            cj = getN(fj)
            cmd = "intersectBed -a %s -b %s -u > ./tmp.bed" % (fi, fj)
            #print(cmd)
            subprocess.call(cmd, shell=True)
            cij = getN("./tmp.bed")
            ji = cij / 1.0 / (ci + cj - cij)
            ds[ni][nj] = ji
    ds = pd.DataFrame(ds)
    ds.to_csv(fnOut, sep="\t", index_label="sample")


def filterPeaks(fin, fout, qcut):
    with open(fout, "w") as fo:
        for line in open(fin):
            line = line.split("\n")[0].split("\t")
            if float(line[8]) > qcut:
                fo.write("\t".join(line) + "\n")


def plotJIClusters(fin, fout):
    ds = pd.read_csv(fin, index_col=0, sep="\t")
    fig, ax = pylab.subplots(figsize=(6, 6))
    g = sns.clustermap(ds.corr(), robust=True, square=True, cmap="vlag")
    #g = sns.clustermap(ds,robust=True,square=True)
    pylab.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
    pylab.setp(g.ax_heatmap.get_xticklabels(), rotation=90)
    #fig.tight_layout()
    fig.set_tight_layout(True)
    pylab.savefig(fout + ".pdf")


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
                jis.append(mat.loc[ns[i], ns[j]])
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
    ax.set_title("Mean JI for reps: %f \n Mean JI for non-reps: %f " %
                 (ds.mean(), np.mean(jis)))
    pylab.savefig(fout + ".pdf")


def main():
    #fs = glob("../../2.peaks/*/*narrowPeak")
    #getJIs(fs, "./peakOverlaps_summary.txt")
    #plotJIClusters("./peakOverlaps_summary.txt","./peakOverlaps")
    evulateReps("./peakOverlaps_summary.txt", "./peakOverlaps_reps")
    #screen for better cutoffs
    """
    rfs = glob("../3.peaks/*.narrowPeak")
    for i in range(2,3):
        nfs = []
        for f in rfs:
            fo = "./tmp/" + f.split("/")[-1]
            nfs.append( fo )
            filterPeaks( f,fo,i)
        getJIs(nfs,"./3.peakOverlaps/3.peakOverlaps_summary_cut%s.txt"%i) 
        evulateReps("./3.peakOverlaps/3.peakOverlaps_summary_cut%s.txt"%i,"./3.peakOverlaps/cut%s"%i) 
    """


if __name__ == "__main__":
    main()
