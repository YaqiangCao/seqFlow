#!/usr/bin/env python2.7
#--coding:utf-8--
"""
scMNase_getFeatureCount.py
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
import numpy as np
import pandas as pd
import HTSeq
from joblib import Parallel, delayed

#this
#from utils import getLogger, callSys, PET
from utils import *

#global settings
#logger
date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
logger = getLogger(fn=os.getcwd() + "/" + date.strip() + "_" +
                   os.path.basename(__file__) + ".log")

global model
model = None


def getModel(modelf):
    model = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    for line in open(modelf):
        line = line.split("\n")[0].split("\t")
        n = "|".join(line)
        iv = HTSeq.GenomicInterval(line[0], int(line[1]), int(line[2]))
        model[iv] = n
    #for iv,v in model.steps():
    #    print(iv,v)
    return model


def getFeatureCount(f,todir="./data",mode="cN"):
    """
    """
    n = f.split("/")[-1].replace(".bedpe.gz", "")
    fout = todir + "/" + n + '_%s.txt'%mode
    if os.path.isfile(fout):
        return
    print(f)
    ss = set()
    for i, line in enumerate(gzip.open(f)):
        if i % 10000 == 0:
            cFlush("%s read from %s" % (i, f))
        line = line.split("\n")[0].split("\t")
        if line[0] != line[3]:
            continue
        s = min(int(line[1]), int(line[4]))
        e = max(int(line[2]), int(line[5]))
        d = e - s 
        m = (s+e)/2
        if mode == "cN" and 140 <= d <= 180:
            #iv = HTSeq.GenomicInterval(line[0], s, e)
            iv = HTSeq.GenomicInterval(line[0], m, m+1)
            for niv, nv in model[iv].steps():
                if nv != set([]):
                    ss.add(nv)
                    #ss.update(nv)
        elif mode == "sP" and d <= 80: 
            #iv = HTSeq.GenomicInterval(line[0], s, e)
            iv = HTSeq.GenomicInterval(line[0], m, m+1)
            for niv, nv in model[iv].steps():
                if nv != set([]):
                    #ss.update(nv)
                    ss.add( nv )
        else:
            continue
    print()
    with open(fout, "w") as fo:
        fo.write("\n".join(list(ss)))
    print(f, "finished")
    logger.info("file:%s,mode:%s,features:%s"%(f,mode,len(ss)))


def summary(pre="pcRNA",suffix="cN",todir="./data"):
    ns = set()
    cs = []
    fs = glob("%s/*%s.txt"%(todir,suffix))
    fs.sort()
    for f in fs:
        n = f.split("/")[-1].split("_%s.txt"%suffix)[0]
        rs = open(f).read().split("\n")
        ns.update(rs)
        cs.append(n)
    ds = np.zeros([len(ns), len(cs)])
    ds = pd.DataFrame(ds, index=ns, columns=cs)
    for f in fs:
        n = f.split("/")[-1].split("_%s.txt"%suffix)[0]
        rs = open(f).read().split("\n")
        ds[n][rs] = 1.0
        print(ds[n][rs])
    ds.to_csv("%s_%s_binary.txt"%(pre,suffix), index_label="pos", sep="\t")


def filterMat(f,cut=20):
    with open(f.replace(".txt","_filter.txt"),"w") as fo:
        for i, line in enumerate(open(f)):
            if i == 0:
                fo.write(line)
                continue
            if i % 1000000 == 0:
                cFlush("%s read from %s"%(i,f))
            line = line.split("\n")[0].split("\t")
            ns = np.array( map(float,line[1:]))
            if np.sum(ns) >= cut:
                fo.write("\t".join(line)+"\n")



def main():
    """
    todir = "./data"
    if not os.path.exists(todir):
        os.mkdir(todir)
    global model
    modelf = "../../../8.SetsForClustering/mm10_whole_genome_wins.bed" 
    model = getModel(modelf)
    fs = glob("../../../5.reduBedpe/*.bedpe.gz")
    ds = pd.read_table("../../../4.bedpe/stat_filter2.txt",index_col=0,sep="\t").index
    fs = [f for f in fs if f.split("/")[-1].split(".bedpe")[0] in ds]
    fs.sort()
    Parallel(n_jobs=10)(delayed(getFeatureCount)(f,todir,"cN") for f in fs)
    Parallel(n_jobs=10)(delayed(getFeatureCount)(f,todir,"sP") for f in fs)
    summary("pcRNA","cN",todir)
    summary("pcRNA","sP",todir)
    """
    for f in glob("*.txt"):
        filterMat(f)


if __name__ == '__main__':
    start_time = datetime.now()
    main()
    elapsed = datetime.now() - start_time
    print "The process is done"
    print "Time used:", elapsed