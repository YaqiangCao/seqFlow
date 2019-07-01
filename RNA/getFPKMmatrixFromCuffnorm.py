#!/usr/bin/env python2.7
#--coding:utf-8--
"""
getFPKMmatrixFromCuffnorm.py
"""

__author__ = "CAO Yaqiang"
__date__ = "2015-03-04"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
from datetime import datetime
import glob, os, time

#3rd library
import matplotlib as mpl
mpl.use("pdf")
mpl.rcParams["pdf.fonttype"] = 42
import pylab
import numpy as np
import pandas as pd
import seaborn as sns
sns.set_style("white")


def getFPKMMatrix(fpkm_table, gene_attr,sample_info ,pre="test"):
    mat = pd.read_table(fpkm_table, index_col=0)
    mat = mat.fillna(0.0)
    #filter all zeros
    nis = []
    for i in mat.index:
        s = mat.loc[i, :]
        if np.sum(s) == 0.0:
            nis += [i]
    mat = mat.drop(nis)
    #assign genes
    genes = pd.read_table(gene_attr, index_col=0)
    ngs = []
    for i in mat.index:
        ni = [i, genes.loc[i, "gene_short_name"], genes.loc[i, "locus"]]
        ni = "|".join(ni)
        ngs += [ni]
    mat.index = ngs
    #assign columns
    ts = {}
    for i,line in enumerate(open(sample_info)):
        if i == 0:
            continue
        line = line.split( "\n" )[ 0 ].split( "\t" )  
        ts[line[0]] = line[1].split("/")[-2]
    ts = pd.Series(ts)
    mat.columns = ts[mat.columns]
    mat.to_csv(pre + ".txt", index_label="gene", sep="\t")


def getFPKMcutoff(f):
    mat = pd.read_table(f, index_col=0)
    mat = mat.fillna(0.0)
    samples = {}
    for c in mat.columns:
        nc = c[-1]
        if nc not in samples:
            samples[nc] = []
        samples[nc] += [c]
    data = []
    for ss, cs in samples.items():
        sub = mat[cs]
        for t in sub.itertuples():
            s = np.array(t[1:])
            ns = s[s > 0.0]
            if len(ns) < len(s) / 2:
                data.extend(list(ns))
        print ss, len(data)
    data = np.array(data)
    data = data[data > 0]
    data = data[~np.isnan(data)]
    data = np.log2(data)
    fig, ax = pylab.subplots()
    ax = sns.kdeplot(data, shade=True)
    ylim = ax.get_ylim()
    cut = np.mean(data) + 3 * np.std(data)
    print(np.mean(data),np.std(data))
    ax.vlines(cut, ylim[0], ylim[1])
    ncut = 2**cut
    #ncut = cut
    ax.text(cut, ylim[1] - 0.1, "cutoff=%s" % ncut)
    ax.set_xlabel("log2(FPKM)")
    ax.set_ylabel("Density")
    #ax.set_xlim( -15,5 )
    pre = f.split("/")[-1].split(".")[0]
    pylab.savefig("%s_FPKMCutoff.pdf" % pre)
    return ncut


def getFilteredMatrix():
    fs = glob.glob("*.txt")
    for f in fs:
        mat = pd.read_table(f, index_col=0)
        #cut = getFPKMcutoff(f)
        cut = 1
        nis = []
        for t in mat.itertuples():
            s = np.array(t[1:]) 
            s = s[s > cut]
            if len(s) == 0:
                nis.append(t[0])
        mat = mat.drop(nis)
        mat = np.log2(mat + 1.0)
        pre = f.split("/")[-1].split(".")[0]
        mat.to_csv(pre + "_filter.txt", sep="\t")


def getGenes(gtf):
    gs = set()
    for line in open(gtf):
        line = line.split( "\n" )[ 0 ].split( "\t" )  
        g = line[8].split(";")[0].split()[1].replace('"','')
        gs.add(g)
    return gs



def main():
    getFPKMMatrix( "cell/genes.fpkm_table","cell/genes.attr_table","cell/samples.table",pre="cell" )
    #getFilteredMatrix()
    #seperateRNAs()


main()
