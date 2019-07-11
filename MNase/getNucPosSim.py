#!/usr/bin/env python
#--coding:utf-8--
"""
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
from utils import *

#global settings
#logger
date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
logger = getLogger(fn=os.getcwd() + "/" + date.strip() + "_" +
                   os.path.basename(__file__) + ".log")

def buildCovModel(readF, dfilter=[80, 140, 180], mapq=1):
    """
    Building Genome Coverage profile for MNase-seq data based on HTSeq.

    Parameters
    ---
    readF: str,bedpe.gz
    dfilter: list, distance to determin conical and particle
    mapq: int, MAPQ cutoff to remove PETs.

    """
    print("building models for %s" % readF)
    n = readF.split('/')[-1].split(".bedpe.gz")[0]
    model = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    reds = set()
    for i, line in enumerate(gzip.open(readF, 'rt')):
        if i % 10000 == 0:
            report = "%s lines genome signal read." % i
            cFlush(report)
        line = line.split("\n")[0].split("\t")
        if len(line) < 7:
            continue
        if line[0] != line[3]:
            continue
        if int(line[7]) < mapq:
            continue
        s = min(int(line[1]), int(line[4]))
        e = max(int(line[2]), int(line[5]))
        d = e - s
        r = (line[0], s, e)
        if r in reds:
            continue
        else:
            reds.add(r)
        m = (s + e) / 2
        #iv = HTSeq.GenomicInterval(line[0], m, m + 1)
        iv = HTSeq.GenomicInterval(line[0], s, e)
        if d >= dfilter[1] and d <= dfilter[2]:
            model[iv] += iv
    return model



def getNucPosSim(fa,fb,fdhs):
    """
    Get the nucleosome position similarities.
    """
    modela = buildCovModel( fa )
    modelb = buildCovModel( fb )
    n = 0
    d = 0
    for line in open(fdhs):
        line = line.split("\n")[0].split("\t")
        iv = HTSeq.GenomicInterval(line[0], int(line[1]),int(line[2]) )
        ra = set()
        for iva,ita in modela[iv].steps():
            if ita != set([]):
                ra.update( ita)
        rb = set()
        for ivb,itb in modelb[iv].steps():
            if itb != set([]):
                rb.update( itb)
        ra = list(ra)
        rb = list(rb)
        #one DHS region may contains multiple nucleosome
        if len(ra) >0 and len(ra) <=3 and len(rb) >0 and len(rb) <=3:
            for a in ra:
                for b in rb:
                    #overlapped nuclosomes
                    if a.overlaps( b ):
                        n += 1
                        nd = abs( (a.start + a.end)/2 - (b.start+b.end)/2 )
                        d += nd
    d = d / 1.0 /n 
    na = fa.split("/")[-1].split(".bedpe.gz")[0]
    nb = fb.split("/")[-1].split(".bedpe.gz")[0]
    print(na,nb,d)
    return na,nb,d


def main():
    fdhs = "Merged_DHSs_mm10.bed"
    ds = pd.read_table("../4.bedpe/stat_filter3.txt",index_col=0,sep="\t").index
    fs = glob("../5.reduBedpe/*.bedpe.gz")
    fs =  [f for f in fs if f.split("/")[-1].split(".bedpe")[0] in ds]
    ps = []
    for i in range(len(fs)):
        for j in range(i+1,len(fs)):
            ps.append( [fs[i],fs[j],fdhs] )
    ds = Parallel(n_jobs=10)(delayed( getNucPosSim )(p[0],p[1],p[2]) for p in ps)
    data = {}
    for d in ds:
        if d[0] not in data:
            data[ d[0] ] = {}
        if d[1] not in data:
            data[ d[1] ] = {}
        data[ d[0] ][ d[1] ] = d[2]
        data[ d[1] ][ d[0] ] = d[2]
    data = pd.DataFrame( data )
    data = data.fillna(0)
    data.to_csv( "NucPosSim.txt",sep="\t")

if __name__ == '__main__':
    start_time = datetime.now()
    main()
    elapsed = datetime.now() - start_time
    print "The process is done"
    print "Time used:", elapsed
