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
import os, time,string
from glob import glob
from datetime import datetime

#3rd library
import numpy as np
import matplotlib as mpl
mpl.use( "pdf" )
import seaborn as sns
mpl.rcParams[ "pdf.fonttype" ] = 42
mpl.rcParams[ "figure.figsize" ] = ( 4,2.75 )
mpl.rcParams[ "figure.dpi" ] = 100
mpl.rcParams[ "savefig.transparent" ] = True
mpl.rcParams[ "savefig.bbox" ] = "tight"
mpl.rcParams[ "font.size" ] = 10.0
mpl.rcParams[ "font.sans-serif" ] = "Arial"
mpl.rcParams[ "savefig.format" ] = "pdf"
import pylab
sns.set_style( "white" )
import numpy as np
import pandas as pd
from numpy.linalg import norm
from scipy.stats import entropy
from joblib import Parallel,delayed



def calcJSD(P,Q):
    _P = P / norm(P, ord=1)
    _Q = Q / norm(Q, ord=1)
    _M = 0.5*(_P + _Q)
    return 1 - 0.5*(entropy(_P,_M) + entropy(_Q,_M))



def defineTs(grp,f):
    ncs = pd.read_table(f,index_col=0).columns 
    ss = pd.Series.from_csv(grp,sep="=")
    ncs = ss.index.intersection(ncs)
    ss = ss[ncs]
    ts = {}
    for s in set(ss.values):
        ns = pd.Series(np.zeros(ss.shape[0]),index=ss.index)
        ns[ss==s] = 1
        ts[s] = ns
    return ts,ss.index



def getJSDs(f,ts,cs):
    mat = pd.read_table(f,index_col=0,sep="\t")
    mat = mat[cs]
    for key in ts.keys():
        if key != "Epithelia":
            continue
        print key
        ss = {}
        for ns in mat.itertuples():
            ss[ns[0]] = calcJSD(ns[1:],ts[key])
        ss = pd.Series(ss)
        ss.to_csv(key+".jsd",sep="\t")


def main():
    f = "../../../11.ELRs_PLRs_Binary_Sets/2.ELRs/ELRs.txt"
    grp = "ELRs.grp"
    ts,cs = defineTs(grp,f)
    getJSDs(f,ts,cs)
    
    
main()
