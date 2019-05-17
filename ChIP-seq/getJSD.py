#!/usr/bin/env python2.7
#--coding:utf-8--
"""
getJSD.py
2016-01-19: caculating jessen-shannon divergence 
2016-07-14: modified JSD caculating and multiple processing
2016-07-18: lncRNA TSS used as comparasion. 
2017-10-30: ELRs, PLRs, non-TE promoters, non-TE enhancers
2019-05-17: slightly modified, to set the pre-defined patterns, make modifications to defineTS
"""

__author__ = "CAO Yaqiang"
__date__ = ""
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import os
from glob import glob

#3rd library
from tqdm import tqdm
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
from numpy.linalg import norm
from scipy.stats import entropy
from joblib import Parallel, delayed


def calcJSD(P, Q):
    _P = P / norm(P, ord=1)
    _Q = Q / norm(Q, ord=1)
    _M = 0.5 * (_P + _Q)
    return 1 - 0.5 * (entropy(_P, _M) + entropy(_Q, _M))


def defineTs(f):
    """
    Define the target pattern. 
    sa: 0 0 0 ... 1 1, which specificy the gain activity in AID, key is -1
    sb: 1 1 1 ... 0 0, which specificy the loss activity in AID, key is 1
    """
    cs = pd.read_table(f, index_col=0).columns
    sa = []
    for c in cs:
        if ("AID" in c) and ("Pos" in c):
            sa.append(0)
        else:
            sa.append(1)
    sa = np.array(sa)
    sb = 1 - sa
    ts = {}
    ts[-1] = sa
    ts[1] = sb
    return ts


def sJSD(ns, ts):
    """
    Caculating JSD for all pre-designed patterns.
    """
    ns = np.array(ns[1:])
    ss = {}
    for key, nt in ts.items():
        d = calcJSD(ns, nt)
        ss[key] = d
    return ss


def evualSpecific(f, ts, fout):
    mat = pd.read_table(f, index_col=0, sep="\t")
    ss = Parallel(n_jobs=1)(delayed(sJSD)(ns, ts)
                            for ns in tqdm(list(mat.itertuples())))
    ds = {}
    for i, n in enumerate(mat.index):
        ds[n] = ss[i]
    ds = pd.DataFrame(ds).T
    ds.to_csv(fout + "_jsd.txt", sep="\t", index_label="peakId")


"""
def plotSpecific(fs):
    fig, ax = pylab.subplots()
    for f in fs:
        label = f.split(".jsd")[0]
        print label
        s = pd.Series.from_csv(f, sep="\t")
        sns.kdeplot(s, label=label, shade=False, cumulative=True)
    ax.legend()
    ax.set_xlabel("Maximal JS Specificity Score")
    ax.set_ylabel("Cumulative Density")
    ax.set_xlim([0.3, 1.0])
    pylab.savefig("jsd.pdf")
"""


def main():
    for f in glob("../10.getNoiseFilter/*.txt"):
        n = f.split("/")[-1].split("_TPM")[0]
        print(n)
        ts = defineTs(f)
        evualSpecific(f, ts, n)
    #plotSpecific(glob("*.jsd"))


main()
