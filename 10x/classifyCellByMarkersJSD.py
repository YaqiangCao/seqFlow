#!/usr/bin/env python
#--coding:utf-8--
"""
classifyCellByMarkersJSD.py

classify the single cells by marker patter JSD.
"""

__author__ = "CAO Yaqiang"
__date__ = ""
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import os, time, string
from glob import glob
from collections import Counter
from datetime import datetime

#3rd library
import numpy as np
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
from tqdm import tqdm
from numpy.linalg import norm
from scipy.stats import entropy
from joblib import Parallel, delayed


def preMat(f):
    """
    Pre-processing the matirx to mean values of each group.
    """
    mat = pd.read_csv(f, index_col=0, sep="\t")
    groups = {}
    for i in mat.index:
        ni = i.split("|")
        ni = ni[0] + "|" + ni[-1]
        if ni not in groups:
            groups[ni] = []
        groups[ni].append(i)
    #summary the mat
    nmat = {}
    for k, v in groups.items():
        smat = mat.loc[v, ]
        smat = smat.mean(axis=0)
        nmat[k] = smat
    nmat = pd.DataFrame(nmat).T
    return nmat


def calcJSD(P, Q):
    _P = P / norm(P, ord=1)
    _Q = Q / norm(Q, ord=1)
    _M = 0.5 * (_P + _Q)
    return 1 - 0.5 * (entropy(_P, _M) + entropy(_Q, _M))


def defineTs(mat):
    """
    Define the target patterns
    """
    ts = {}
    ss = {i: i for i in mat.index}
    ss = pd.Series(ss)
    for k in set(ss.values):
        ins = ss[ss == k].index
        ns = pd.Series(np.zeros(ss.shape[0]), index=ss.index)
        ns[ins] = 1
        ts[k] = ns
    return ts


def getJSD(mat, ts, pre):
    """
    Get the JSD for each cell.
    """
    rs = {}
    for c in tqdm(mat.columns):
        sa = mat[c]
        ss = {}
        for key in ts.keys():
            sb = ts[key]
            ss[key] = calcJSD(sb, sa)
        ss = pd.Series(ss)
        rs[c] = ss
    rs = pd.DataFrame(rs).T
    rs.to_csv(pre + "_jsd.txt", sep="\t")


def getCellTypes(f):
    mat = pd.read_csv(f,index_col=0,sep="\t")
    rs = []
    for t in mat.itertuples():
        p = np.argmax( t[1:])
        rs.append( mat.columns[p] )
    rs = pd.Series( Counter(rs))
    return rs


def main():
    """
    f = "../4.showMarker/GR_3_marker.txt"
    mat = preMat(f)
    ts = defineTs(mat)
    getJSD(mat, ts, "GR_3")
    """
    ds = {}
    for f in sorted(glob("*.txt")):
        n = f.split("_jsd")[0]
        rs = getCellTypes(f)
        ds[n] = rs
    ds = pd.DataFrame(ds)
    ds = ds.fillna(0)
    ds.to_csv("test.txt",sep="\t")


main()
