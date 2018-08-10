#!/usr/bin/env python2.7
#--coding:utf-8--
"""
selJSD.py
2016-01-19: caculating jessen-shannon divergence 
2016-07-14: modified JSD caculating and multiple processing
2016-07-18: lncRNA TSS used as comparasion. 
2017-10-30: ELRs, PLRs, non-TE promoters, non-TE enhancers
"""

__author__ = "CAO Yaqiang"
__date__ = ""
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import os, time, string
from glob import glob
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
from numpy.linalg import norm
from scipy.stats import entropy
from joblib import Parallel, delayed


def calcJSD(P, Q):
    _P = P / norm(P, ord=1)
    _Q = Q / norm(Q, ord=1)
    _M = 0.5 * (_P + _Q)
    return 1 - 0.5 * (entropy(_P, _M) + entropy(_Q, _M))


def defineTs(grp, f):
    ncs = pd.read_table(f, index_col=0).columns
    ss = pd.Series.from_csv(grp, sep="=")
    ncs = ss.index.intersection(ncs)
    ss = ss[ncs]
    ts = {}
    for s in set(ss.values):
        ns = pd.Series(np.zeros(ss.shape[0]), index=ss.index)
        ns[ss == s] = 1
        ts[s] = ns
    return ts, ss.index


def sJSD(ns, ts):
    ns = np.array(ns[1:])
    tmp = 0
    for nt in ts.values():
        d = calcJSD(ns, nt)
        if d > tmp:
            tmp = d
    return tmp


def evualSpecific(f, cs, ts, fout):
    print f
    mat = pd.read_table(f, index_col=0, sep="\t")
    mat = mat[cs]
    print mat.index[:10]
    print mat.shape
    ns = [t[0] for t in mat.itertuples() if np.sum(t[1:]) <= 1]
    mat = mat.drop(ns)
    print mat.shape
    ss = Parallel(n_jobs=1)(delayed(sJSD)(ns, ts) for ns in mat.itertuples())
    ss = pd.Series(ss, index=mat.index)
    ss.to_csv(fout + ".jsd", sep="\t")


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


def main():
    """
    grp = "ELRs.grp"
    ts,cs = defineTs(grp,"../11.ELRs_PLRs_Binary_Sets/2.ELRs/ELRs.txt")
    
    evualSpecific("../11.ELRs_PLRs_Binary_Sets/2.ELRs/ELRs.txt",cs,ts,"ELRs")
    evualSpecific("../11.ELRs_PLRs_Binary_Sets/1.PLRs/PLRs.txt",cs,ts,"PLRs")
    evualSpecific("../11.ELRs_PLRs_Binary_Sets/5.Non-TEs_TypicalPromoters/2.BinarySet/promoters.mat",cs,ts,"noTEPromoters")
    evualSpecific("../11.ELRs_PLRs_Binary_Sets/6.Non-TEs_TypicalEnhancers/2.BinarySet/enhancers.mat",cs,ts,"noTEEnhancers")
    evualSpecific("../11.ELRs_PLRs_Binary_Sets/3.TypicalPromoters/2.BinarySet/promoters.mat",cs,ts,"typicalPromoters")
    evualSpecific("../11.ELRs_PLRs_Binary_Sets/4.TypicalEnhancers/2.BinarySet/enhancers.mat",cs,ts,"typicalEnhancers")
    """
    plotSpecific(glob("*.jsd"))


main()
