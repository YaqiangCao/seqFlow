#!/usr/bin/env python2.7
#--coding:utf-8--
"""
carrykmeans.py
2016-01-15: specific analysis for binary matrix  
"""

__author__ = "CAO Yaqiang"
__date__ = "2016-01-15"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import glob, os, time, string
from datetime import datetime

#3rd library
import matplotlib as mpl
mpl.use("pdf")
import seaborn as sns
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["figure.figsize"] = (4, 2.75)
mpl.rcParams["figure.dpi"] = 100
mpl.rcParams["savefig.transparent"] = True
mpl.rcParams["savefig.bbox"] = "tight"
mpl.rcParams["font.size"] = 10.0
mpl.rcParams["font.sans-serif"] = "Arial"
mpl.rcParams["savefig.format"] = "pdf"
import pylab
sns.set_style("white")
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn import metrics
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

#plot setting


def txt2cdt(f):
    mat = pd.read_table(f, index_col=0)
    genes, mat, columns = list(mat.index), mat.values, list(mat.columns)
    #genes = [ g.split( "|" )[ 0 ] for g in genes ]
    genes = map(str, genes)
    nf = f.replace(".txt", ".cdt")
    with open(nf, "w") as f:
        column_header = string.join(['UNIQID', 'NAME', 'GWEIGHT'] + columns,
                                    '\t') + '\n'
        f.write(column_header)
        eweight = string.join(
            ['EWEIGHT', '', ''] + ['1'] * len(columns),
            '\t') + '\n'  ### format column-flat-clusters for export
        f.write(eweight)
        for i, g in enumerate(genes):
            line = string.join([g] * 2 + ['1'] + map(str, mat[i]), '\t') + '\n'
            f.write(line)


def getcs(grp="ELRs.grp"):
    cs = pd.Series.from_csv(grp, sep="=")
    return cs.index


def evualk(f, k=90):
    cs = getcs()
    mat = pd.read_table(f, index_col=0)
    mat = mat[cs]
    ns = [t[0] for t in mat.itertuples() if np.sum(t[1:]) <= 1]
    mat = mat.drop(ns)
    pca = PCA(n_components=k).fit(mat)
    rs = pca.explained_variance_ratio_ * 100
    #rs = np.cumsum(rs)
    fig, ax = pylab.subplots()
    ks = np.arange(1, k + 1)
    ax.plot(ks, rs, marker="o")
    ax.set_xlabel("PCs")
    ax.set_ylabel("Explained Variance Ratio")
    pylab.savefig("pcs.pdf")


def ck(f, pre, k):
    cs = getcs()
    mat = pd.read_table(f, index_col=0)
    mat = mat[cs]
    pca = PCA(n_components=k).fit(mat)
    ns = [t[0] for t in mat.itertuples() if np.sum(t[1:]) <= 1]
    mat = mat.drop(ns)
    c = KMeans(
        n_clusters=k,
        init=pca.components_,
        n_init=1,
        n_jobs=20,
        random_state=123).fit(mat)
    c.fit(mat)
    ds = pd.Series(c.labels_, index=mat.index)
    cs = []
    for t in sorted(set(c.labels_)):
        s = ds[ds == t]
        cs.extend(list(s.index))
    mat = mat.loc[cs, :]
    fout = "%s.txt" % pre
    mat.to_csv(fout, sep="\t", index_label="rep")
    ds.sort()
    ds.to_csv(pre + ".kgg", sep="\t", index_label="rep")
    txt2cdt(fout)


def main():
    f = "../../11.ELRs_PLRs_Binary_Sets/2.ELRs/ELRs.txt"
    evualk(f, k=30)
    ck(f, "ELRs_k12", 12)
    ck(f, "ELRs_k13", 13)
    ck(f, "ELRs_k14", 14)
    ck(f, "ELRs_k15", 15)


main()
