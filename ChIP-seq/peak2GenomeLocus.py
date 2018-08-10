#!/usr/bin/env python2.7
#--coding:utf-8 --
"""
"""

__author__ = "CAO Yaqiang"
__date__ = "2015-03-25"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import glob, os, time
from datetime import datetime

#3rd library
import matplotlib as mpl
mpl.use("pdf")
mpl.rcParams["pdf.fonttype"] = 42
import pandas as pd
import numpy as np
import pylab
#import seaborn as sns
import prettyplotlib as ppl
import brewer2mpl
colors = brewer2mpl.get_map('Set3', 'qualitative', 8).mpl_colors


def read_stat(pattern="", pre="test"):
    sts = glob.glob("../1.*/*%s*.stat" % pattern)
    data = {}
    for s in sts:
        sample = s.split("/")[-1].split("_peaks")[0]
        #print sample,factor
        if sample not in data:
            data[sample] = {
                "3UTR": 0,
                "5UTR": 0,
                "Exon": 0,
                "Intron": 0,
                "Intergenic": 0,
                "Promoter": 0,
                "TTS": 0,
                "Others": 0
            }
        #only collect 3UTR 5UTR Exon Intron Intergenic, TTS, Promoter
        #peaks means the peaks number, bps means the genome occupancy
        ns = open(s).read().split("\nAnnotation")[0].split("\n")[1:]
        for n in ns:
            n = n.split("\t")
            if n[0] in data[sample]:
                data[sample][n[0]] = float(n[1])
            else:
                data[sample]["Others"] += float(n[1])
    #genomic size for each features
    ref = {
        "3UTR": 0,
        "5UTR": 0,
        "Exon": 0,
        "Intron": 0,
        "Intergenic": 0,
        "Promoter": 0,
        "TTS": 0,
        "Others": 0
    }
    ns = open(s).read().split("\nAnnotation")[0].split("\n")[1:]
    for n in ns:
        n = n.split("\t")
        if n[0] in ref:
            ref[n[0]] = float(n[2])
        else:
            ref["Others"] += float(n[2])
    data["ref"] = ref
    data = pd.DataFrame(data)
    data.to_csv(pre + "_features_stats.txt", sep="\t", index_label="feature")


def plotRatio(f="GenomicDis.txt"):
    mat = pd.read_table(f, index_col=0)
    for c in mat.columns:
        mat[c] = mat[c] / np.sum(mat[c])
    mat = mat * 100
    nis = []
    for i in mat.index:
        s = mat.loc[i]
        s = s[s > 5]
        if len(s) == 0:
            nis.append(i)
    mat = mat.drop(nis)
    f, ax = pylab.subplots(figsize=(4, 6))
    x = np.arange(len(mat.columns))
    y = np.zeros(len(mat.columns))
    for i, t in enumerate(mat.index):
        ax.bar(x,
               mat.loc[t, :],
               bottom=y,
               label=t,
               alpha=0.9,
               color=colors[i],
               edgecolor="white")
        y += mat.loc[t, :]
    ax.set_xticks(x)
    ax.set_xticklabels(list(mat.columns), rotation=45)
    ax.set_ylabel("Percentage %")
    leg = ax.legend(loc="upper left", fancybox=True, bbox_to_anchor=(1, 1))
    pylab.savefig("GenomicDis.pdf")


def main():
    #read_stat( pattern="SMAD2",pre="SMAD2" )
    plotRatio(f="SMAD2_features_stats.txt")


if __name__ == "__main__":
    main()
