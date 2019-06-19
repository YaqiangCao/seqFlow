#!/usr/bin/env python2.7
#--coding:utf-8--
"""
callSplicing_rMAST.py
"""

__author__ = "CAO Yaqiang"
__date__ = ""
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import os, time
from glob import glob
from datetime import datetime

#3rd library
import numpy as np
import pandas as pd
import seaborn as sns
from joblib import Parallel, delayed
#my own
from Biolibs.rel.General.logger import getlogger

#global settings
#logger
date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
logger = getlogger(fn=os.getcwd() + "/" + date.strip() + "_" +
                   os.path.basename(__file__) + ".log")

#global
GTF = "/picb/molsysbio2/caoyaqiang/2.WWX_RNA-seq/1.Reference/2.Genecode/2.Parsed/genecode.v27.gtf"
rmast = "/picb/molsysbio2/caoyaqiang/2.WWX_RNA-seq/0.Tools/rMATS.4.0.1/rMATS-turbo-Linux-UCS4/rmats.py"


def call_sys(cmds):
    """
    Call systematic commands without return.
    """
    for c in cmds:
        logger.info(c)
        try:
            os.system(c)
        except:
            logger.error(c)


def preDs(fs):
    ds = {}
    for f in fs:
        n = f.split("/")[-2]
        nn = n[:-1]
        if nn not in ds:
            ds[nn] = {}
        ds[nn][n[-1]] = f
    nds = []
    for key in ds.keys():
        if len(ds[key]) != 2:
            continue
        t = (key, ds[key]["L"], ds[key]["T"])
        nds.append(t)
    return nds


def callSplicing(t):
    out = t[0]
    if os.path.exists(out):
        print(out + "exists, return")
        return
    ta = t[1].split("/")[-1]
    tb = t[2].split("/")[-1]
    c1 = "echo %s > %s.txt" % (t[1], ta)
    c2 = "echo %s > %s.txt" % (t[2], tb)
    c3 = "python2.7 {rmast} --gtf {gtf} --b1 {b1} --b2 {b2} --od {out} -t paired --readLength 150 --nthread 10 --tstat 10".format(
        rmast=rmast, gtf=GTF, b1=ta + ".txt", b2=tb + ".txt", out=out)
    c4 = "rm %s.txt %s.txt" % (ta, tb)
    call_sys([c1, c2, c3, c4])


def summarySplicing(r_cut=20):
    """
    using files with suffix as JCEC,IncLevel could be used,also filtering reads numbers
    """
    terms = ["SE", "RI", "MXE", "A3SS", "A5SS"]
    suffix = ".MATS.JCEC.txt"
    for term in terms:
        fout = term + "_summary.txt"
        if os.path.isfile(fout):
            continue
        fs = glob("*/%s" % (term + suffix))
        ds = {}
        for f in fs:
            print(f)
            n = f.split("/")[-2]
            na = n + "L"
            nb = n + "T"
            sa = {}
            sb = {}
            mat = pd.read_table(f, index_col=0)
            for t in mat.itertuples():
                gid = "|".join(map(str, t[1:11]))
                nt = np.array(t[12:16])
                nt = nt[nt > r_cut]
                if len(nt) < 1:
                    continue
                sa[gid] = t[20]
                sb[gid] = t[21]
            ds[na] = sa
            ds[nb] = sb
        ds = pd.DataFrame(ds)
        ds.to_csv(fout, sep="\t", index_label="id")


def filterSummary(f, n_cut=50):
    ds = pd.read_table(f, index_col=0)
    print(ds.shape)
    ns = []
    for t in ds.itertuples():
        nt = pd.Series(t)
        nt = pd.isnull(nt)
        nt = nt[nt == False]
        if len(nt) < n_cut:
            ns.append(t[0])
    ds = ds.drop(ns)
    print(ds.shape)
    ds.to_csv(f, sep="\t")


def main():
    fs = glob("../2.Mapping/*/*.bam")
    ds = preDs(fs)
    #Parallel(n_jobs=20)(delayed(callSplicing)(t) for t in ds)
    summarySplicing()
    map(filterSummary, glob("*.txt"))


main()
