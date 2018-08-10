#!/usr/bin/env python2.7
#--coding:utf-8--
"""
getFPM.py
2015-11-30: modified pd.DataFrame itertuples
2015-12-18: modified for roadmap.
"""

__author__ = "CAO Yaqiang"
__date__ = ""
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import os, glob, gzip, sys, random
from datetime import datetime
from collections import Counter

#plotting settings

#3rd library
import HTSeq
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

#global setting


def commandFlush(r):
    """
    One line flush.
    """
    sys.stdout.write("\r%s" % r)
    sys.stdout.flush()


def bg2GModel(bg):
    """
    BedGraph format, gzip or not into HTSeq.GenomicArray 
    """
    if bg.endswith(".gz"):
        f = gzip.open(bg, "rb")
    else:
        f = open(bg)
    print datetime.now(), "Start building model for %s" % bg
    model = HTSeq.GenomicArray("auto", stranded=False)
    for i, line in enumerate(f):
        if i % 10000 == 0:
            report = "%s lines genome signal read." % i
            commandFlush(report)
        line = line.split("\n")[0].split("\t")
        if len(line) < 3:
            continue
        chrom = line[0]
        s = int(line[1])
        e = int(line[2])
        iv = HTSeq.GenomicInterval(chrom, s, e)
        model[iv] = float(line[3])
    print
    print datetime.now(), "Model built for %s" % bg
    #return genomic coverage model, chromosomes, reads count and read length
    return model


def readRep(repf):
    """
    Read repeats and add iv.
    """
    print datetime.now(), "Start to read all repeats from %s" % repf
    mat = pd.read_table(repf, index_col=0)
    ds = {}
    i = 0
    #for r,t in mat.iterrows(  ):
    for t in mat.itertuples():
        if i % 10000 == 0:
            report = "%s repeats read and processed" % i
            commandFlush(report)
        i += 1
        #iv = HTSeq.GenomicInterval( t[ "chr" ],t[ "s" ],t[ "e" ] )
        iv = HTSeq.GenomicInterval(t[1], t[2], t[3])
        ds[t[0]] = iv
    ds = pd.Series(ds)
    print
    return ds


def getCoverage(model, iv):
    """
    Get reads count from model in defined region.
    """
    c = 0
    for ivb, value in model[iv].steps():
        c += ivb.length * value
    c = c / iv.length
    return c


def getFPM(sbg, repF):
    reps = readRep(repF)
    pre = sbg.split("/")[-1].split("_")[-2]
    print datetime.now(), "starting getting FPM for %s" % pre
    sModel = bg2GModel(sbg)
    ds = {}
    for i in reps.index:
        ds[i] = getCoverage(sModel, reps[i])
    ds = pd.Series(ds)
    #ds = ds[ ds>0 ]
    ds = ds / ds.sum() * 10**6
    print datetime.now(), "getting FPM for %s finished" % pre
    return pre, ds


def getAllFPM(bgs, repF, pre):
    fout = "%s_FPM.txt" % pre
    if os.path.exists(fout):
        print "%s has been generated, return" % fout
        return
    data = Parallel(n_jobs=len(bgs))(delayed(getFPM)(bg, repF) for bg in bgs)
    mat = {}
    for d in data:
        mat[d[0]] = d[1]
    mat = pd.DataFrame(mat)
    ss = []
    for t in mat.itertuples():
        if np.sum(t[1:]) == 0:
            ss.append(t[0])
    mat = mat.drop(ss)
    mat.to_csv(fout, sep="\t", index_label="rep")


def main():
    repF = "/picb/molsysbio/usr/caoyaqiang/1.Projects/17.Alu/1.GenomeReference/1.hg38/6.Repeats_Annotation/1.GenomicLocation/hg38Rep.txt"
    root = "../../4.NormalizationFactor/1.Consolidated"
    rs = glob.glob(root + "/*.txt")
    random.shuffle(rs)
    for r in rs:
        cell = r.split("/")[-1].replace("_estR.txt", "")
        bgs = glob.glob("%s/%s*.bg.gz" % (root, cell))
        getAllFPM(bgs, repF, cell)


if __name__ == '__main__':
    start_time = datetime.now()
    main()
    elapsed = datetime.now() - start_time
    fn = os.path.basename(__file__)
    print datetime.now(), "The process is done for %s,time used:%s" % (fn,
                                                                       elapsed)
