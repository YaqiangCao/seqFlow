#!/usr/bin/env python2.7
#--coding:utf-8--
"""
getEnrichedTEs.py
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
from scipy.stats import poisson

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
    return i, model


def getCount(t, model, iv):
    ss = set()
    for niv, nv in model[iv].steps():
        if nv != set([]):
            ss.update(nv)
    c = len(ss)
    rpkm = c / 1.0 / iv.length / t * 10**9
    return c, rpkm


def countTEs(f, repf, fout, psedo=1, ext=5):
    """
    Count reads located in TEs and get their enrichment.
    """
    t, model = getCov(f)
    reps = pd.read_table(repf, index_col=0, sep="\t")
    ds = {}
    for rep in tqdm(list(reps.itertuples())):
        rid = rep[0]
        iv = HTSeq.GenomicInterval(rep[1], rep[2], rep[3])
        c, rpkm = getCount(t, model, iv)
        if c == 0:
            continue
        upiv = HTSeq.GenomicInterval(rep[1], rep[2] - iv.length * ext, rep[2])
        upc, uprpkm = getCount(t, model, upiv)
        downiv = HTSeq.GenomicInterval(rep[1], rep[3],
                                       rep[3] + iv.length * ext)
        downc, downrpkm = getCount(t, model, downiv)
        if upc+downc >0:
            es = c / 1.0 / (upc+downc) * 2 * ext
            p = max([1e-300, poisson.sf( c, (upc+downc) / 2.0/ext )])
        else:
            es = c /1.0 / psedo 
            p - 1e-300
        ds[rid] = {
            "length": iv.length,
            "count": c,
            "RPKM": rpkm,
            "up_count_ext%s"%ext: upc,
            "up_RPKM_ext%s"%ext: uprpkm,
            "down_count_ext%s"%ext: downc,
            "down_RPKM_ext%s"%ext: downrpkm,
            "ES": es,
            "poisson_p-value":p,
            #"ES": rpkm / 1.0 /
            #(uprpkm + downrpkm + psedo) * 2,  #psedo count to avoid divid zero
        }
    ds = pd.DataFrame(ds).T
    ds.to_csv(fout + ".txt", sep="\t")


def filterTEs(escut=2.0,pcut=1e-3,lencut=1000):
    ds = set()
    for f in glob("*.txt"):
        mat = pd.read_table(f,index_col=0)
        es = mat["ES"]
        es = es[es>=escut]
        length = mat.loc[es.index,"length"]
        length = length[length>=lencut]
        p = mat.loc[length.index,"poisson_p-value"]
        p = p[p<pcut]
        #ds.update( length.index)
        with open(f.replace(".txt",".bed"),"w") as fo:
            for d in length.index:
                nd = d.split("|")
                line = [nd[0],nd[1],nd[2],d]
                fo.write( "\t".join(line)+"\n")
    #print(ds)
    """
    with open("enrichedTEs.bed","w") as fo:
        for d in ds:
            nd = d.split("|")
            line = [nd[0],nd[1],nd[2],d]
            fo.write( "\t".join(line)+"\n")
    """
            

def main():
    countTEs("test/test.bedpe.gz","test/testRep.txt","test_p") 
    filterTEs()
    #repf = "/home/caoy7/caoy7/Projects/0.Reference/2.mm10/5.repeats/mm10Reps.txt"
    #fs = glob("../3.sepNcSpBedpe/*cN.bedpe.gz")
    #Parallel(n_jobs=len(fs))(delayed(countTEs)(f,repf,f.split("/")[-1].split(".bedpe")[0]+"_reps") for f in fs)
    #filterTEs()

if __name__ == '__main__':
    start_time = datetime.now()
    main()
    elapsed = datetime.now() - start_time
    print "The process is done"
    print "Time used:", elapsed
