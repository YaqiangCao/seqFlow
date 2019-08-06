#!/usr/bin/env python
#--coding:utf-8--
"""
filterCellsByTssProfile.py
"""

__author__ = "CAO Yaqiang"
__date__ = "2019-05-29"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import os, time, commands, gzip
from glob import glob
from datetime import datetime
from collections import Counter

#3rd library
import HTSeq
import numpy as np
import brewer2mpl
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed
#global settings
import matplotlib as mpl
mpl.use("pdf")
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["figure.figsize"] = (4*0.8, 2.75*0.8)
mpl.rcParams["font.size"] = 10.0
from mpl_toolkits.axes_grid.inset_locator import inset_axes
colors = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors
import pylab
import seaborn as sns
sns.set_style("white")


#this
from utils import *

#global settings
#logger
date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
logger = getLogger(fn=os.getcwd() + "/" + date.strip() + "_" +
                   os.path.basename(__file__) + ".log")

def binData(s,start=-1000,end=1000,bin=10):
    s = s.reshape((len(s) / bin, bin))
    s = np.mean(s, axis=1)
    s = pd.Series(s, index=np.arange(start, end, bin)) 
    return s

def checkPattern(fcn,fsp,rcn,rsp):
    ns = rcn.index
    ns = ns[ns>=0]
    tcn = binData( pd.Series.from_csv( fcn,sep="\t").values )
    tsp = binData( pd.Series.from_csv( fsp,sep="\t").values )
    cncorr = tcn.corr( rcn ) 
    spcorr = tcn.corr( rsp )
    corr = tcn.corr(tsp)
    return cncorr,spcorr,corr


def getAllPatternCorr():
    cn = pd.Series.from_csv("cN.txt",sep="\t")
    sp = pd.Series.from_csv("sP.txt",sep="\t")
    ns = cn.index
    ns = ns[ ns >=-1000]
    ns = ns[ ns < 1000]
    cn = cn[ ns ]
    sp = sp[ ns ]
    
    nfs = {}
    for f in glob("../data/*.txt"):
        t = f.split("_")[-1].split(".txt")[0]
        nf ="_".join( f.split("/")[-1].split("_")[:-1] )
        if nf not in nfs:
            nfs[nf] = {}
        nfs[nf][t] =f 
    ds = {}
    for n in tqdm(nfs.keys()):
        cncorr,spcorr,corr = checkPattern( nfs[n]["cn"],nfs[n]["sp"],cn,sp)
        ds[ n ] = {
            "cN_RefcN_PCC": cncorr,
            "cN_RefsP_PCC": spcorr,
            "cN_sP_PCC": corr
            }
    ds = pd.DataFrame(ds).T
    ds.to_csv("cN_cN_sP_corr.txt",sep="\t")


def plotPatternCorr():
    mat = pd.read_table("cN_cN_sP_corr.txt",index_col=0,sep="\t")
    cs = mat.index
    ncs = ["_".join(c.split("_")[:-1]) for c in cs]
    nc = list(set(ncs))
    nc = { c:i for i,c in enumerate(nc) }
    ncs = np.array([nc[c] for c in ncs])

    fig, ax = pylab.subplots()
    x = mat["cN_RefcN_PCC"]
    y = mat["cN_RefsP_PCC"]
    for label, c in nc.items():
        ns = np.where(ncs==c)[0]
        ax.scatter(x[ns], y[ns], color=colors[c], label=label,s=2)
        ax.scatter(1, 1, color=colors[c], label=label, s=0)
    ax.set_xlabel( "cN corr merged all cN")
    ax.set_ylabel( "cN corr merged all sP")
    pylab.tight_layout()
    pylab.savefig("cN_cN_sP_corr.pdf")

def main():
    #getAllPatternCorr()
    plotPatternCorr()

if __name__ == '__main__':
    start_time = datetime.now()
    main()
    elapsed = datetime.now() - start_time
    print "The process is done"
    print "Time used:", elapsed
