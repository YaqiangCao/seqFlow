#!/usr/bin/env python2.7
#--coding:utf-8--
"""
getFRiP.py
Get the ratio of reads that in peaks.
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


def getCov(f,paired=True):
    logger.info("Building coverage model for %s, paired=%s"%(f,paired))
    model = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    i = None
    uniqs = set()
    for i, line in enumerate(gzip.open(f)):
        if i % 10000 == 0:
            cFlush("%s read from %s" % (i, f))
        line = line.split("\n")[0].split("\t")
        if paired and line[0] != line[3]:
                continue
        if paired:
            s = min(int(line[1]), int(line[4]))
            e = max(int(line[2]), int(line[5]))
        else:
            s = int(line[1])
            e = int(line[2])
        r = (line[0],s,e)
        if r not in uniqs:
            iv = HTSeq.GenomicInterval(line[0], s, e)
            model[iv] += str(i)
            uniqs.add( r )
    if i is None:
        logger.error("ERROR! No read in %s."%f)
        return 0, None
    logger.info("%s read from %s, unique %s"%(i,f,len(uniqs)))
    return len(uniqs), model


def getCount(t, model, iv):
    ss = set()
    for niv, nv in model[iv].steps():
        if nv != set([]):
            ss.update(nv)
    c = len(ss)
    rpkm = c / 1.0 / iv.length / t * 10**9
    return c, rpkm


def countFeatures(f, featuref,paired=True):
    """
    Count reads enrichment at the TSS regions.
    """
    #logger.info("Building coverage for %s"%f)
    t, model = getCov(f,paired=paired)
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
    logger.info("FLAG!\t %s:: total:%s,inPeaks:%s,inPeaksRatio:%s"%(n,t,r,r/1.0/t))
    return n, t, r, r/1.0/t
        

def main():
    ds = []
    for f in glob("./*/*.narrowPeak"):
        pre = f.split("/")[-2]
        #bed = "../3.beds/" + pre + ".bed.gz"
        bed = "../3.bedpe/" + pre + ".bedpe.gz"
        if not os.path.isfile(bed):
            continue
        ds.append([bed, f])
    cpu = min(30,len(ds))
    rs = Parallel(n_jobs=cpu)(delayed(countFeatures)(d[0],d[1],paired=True) for d in ds)
    #rs = Parallel(n_jobs=1)(delayed(countFeatures)(d[0],d[1],paired=True) for d in ds)
    ds = {}
    for r in rs:
        ds[ r[0] ] = {
            "uniqueReads": r[1],
            "readsInFeature": r[2],
            "readsRatioInFeature":r[3]
            }
    ds = pd.DataFrame(ds).T
    ds.to_csv("FRiP.txt",sep="\t",index_label="sample")


if __name__ == '__main__':
    start_time = datetime.now()
    main()
    elapsed = datetime.now() - start_time
    print "The process is done"
    print "Time used:", elapsed
