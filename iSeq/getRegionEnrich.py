#!/usr/bin/env python2.7
#--coding:utf-8--
"""
getTssSigEnrich.py
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
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed

#this
#from utils import getLogger, callSys, PET
from utils import *

#global settings
#logger
date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
logger = getLogger(fn=os.getcwd() + "/" + date.strip() + "_" +
                   os.path.basename(__file__) + ".log")


def getCov(f):
    model = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    i = None
    for i, line in enumerate(gzip.open(f)):
        if i % 10000 == 0:
            cFlush("%s read from %s" % (i, f))
        line = line.split("\n")[0].split("\t")
        if line[0] != line[3]:
            continue
        s = min(int(line[1]), int(line[4]))
        e = max(int(line[2]), int(line[5]))
        m = (s + e) / 2
        iv = HTSeq.GenomicInterval(line[0], m, m + 1)
        model[iv] += str(i)
    if i is None:
        return 0, None
    return i, model


def getCount(t, model, iv):
    ss = set()
    for niv, nv in model[iv].steps():
        if nv != set([]):
            ss.update(nv)
    c = len(ss)
    rpkm = c / 1.0 / iv.length / t * 10**9
    return c, rpkm


def countFeatures(f, featuref):
    """
    Count reads enrichment at the TSS regions.
    """
    logger.info("Building coverage for %s"%f)
    t, model = getCov(f)
    if t == 0:
        return None
    logger.info("%s reads from %s"%(t,f))
    ds = {}
    logger.info("Caculating enriched reads at target regions of %s"%featuref)
    r = 0
    for line in tqdm(open(featuref).read().split("\n")):
        line = line.split("\n")[0].split("\t")
        if len(line) < 3:
            continue
        s = int(line[1]) 
        e = int(line[2])
        iv = HTSeq.GenomicInterval(line[0],s,e)
        c, rpkm = getCount(t, model, iv) 
        r += c
    n = f.split("/")[-1].split(".bedpe")[0]
    print(n,t,r)
    return n, t, r, r/1.0/t
        

def main():
    """
    href="/home/caoy7/caoy7/Projects/4.indexing/1.20190607_KZ1822/0.Pub/1.ENCODE_mouseENCODE/GM12878_DHSs.bed"
    fs = glob("../1.human/*.bedpe.gz")
    rs =Parallel(n_jobs=20)(delayed( countFeatures )( f,href ) for f in fs)
    ds = {}
    for r in rs:
        if r is None:
            continue
        ds[r[0]] = {
            "total": r[1],
            "inDHS": r[2],
            "inDHSRatio": r[3],
            }
    ds = pd.DataFrame(ds).T
    ds.to_csv("human_sc_stat.txt",sep="\t")
    href="/home/caoy7/caoy7/Projects/4.indexing/1.20190607_KZ1822/0.Pub/1.ENCODE_mouseENCODE/NIH3T3_DHSs.bed"
    fs = glob("../2.mouse/*.bedpe.gz")
    rs =Parallel(n_jobs=20)(delayed( countFeatures )( f,href ) for f in fs)
    ds = {}
    for r in rs:
        if r is None:
            continue
        ds[r[0]] = {
            "total": r[1],
            "inDHS": r[2],
            "inDHSRatio": r[3],
            }
    ds = pd.DataFrame(ds).T
    ds.to_csv("mouse_sc_stat.txt",sep="\t")
    """
    math = pd.read_table("human_sc_stat.txt",sep="\t",index_col=0)
    matm = pd.read_table("mouse_sc_stat.txt",sep="\t",index_col=0)
    ds = { 
        "humanTotal": math["total"],
            "humanInDHS":math["inDHS"],
            "humanInDHSRatio": math["inDHSRatio"],
        "mouseTotal": matm["total"],
            "mouseInDHS":matm["inDHS"],
            "mouseInDHSRatio": matm["inDHSRatio"],
        }
    ds = pd.DataFrame(ds)
    ds = ds.fillna(0)
    ds.to_csv("human_mouse_in_DHS.txt",sep="\t")




if __name__ == '__main__':
    start_time = datetime.now()
    main()
    elapsed = datetime.now() - start_time
    print "The process is done"
    print "Time used:", elapsed
