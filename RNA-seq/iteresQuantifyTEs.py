#!/usr/bin/env python2.7
#--coding:utf-8--
"""
run_iteres.py
"""

__author__ = "CAO Yaqiang"
__date__ = "2014-11-16"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import glob, os, time, commands, shutil
from datetime import datetime

#3rd library
import pysam, HTSeq
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import pylab
import brewer2mpl

#my own
from utils import getlogger, call_sys

#global settings
#logger
date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
logger = getlogger(fn=os.getcwd() + "/" + date.strip() + "_" +
                   os.path.basename(__file__) + ".log")
#data
REPEAT_GTF = "/picb/molsysbio/usr/caoyaqiang/1.Projects/17.Alu/1.GenomeReference/1.hg38/6.Repeats_Annotation/hg38_sorted_repeatmasker.gtf"
CHROM_SIZE = "/home/caoyaqiang/caoyaqiang_han/1.Projects/17.Alu/1.GenomeReference/1.hg38/1.hg38_Sequence/hg38.chrom.sizes"
RMSK = "/picb/molsysbio/usr/caoyaqiang/1.Projects/17.Alu/1.GenomeReference/1.hg38/6.Repeats_Annotation/rmsk.txt"
SUBFAMILY_SIZE = "/picb/molsysbio/usr/caoyaqiang/1.Projects/17.Alu/1.GenomeReference/1.hg38/6.Repeats_Annotation/subfam.size"
REPEATS_LOCUS = "/picb/molsysbio/usr/caoyaqiang/1.Projects/17.Alu/1.GenomeReference/1.hg38/6.Repeats_Annotation/1.GenomicLocation/hg38_rep_locus.txt"


def iteres_stat_filter(bam,
                       stat_root="0.iteres/1.stat/",
                       filter_root="0.iteres/2.filter/"):
    sample = os.path.splitext(os.path.split(bam)[1])[0]
    if not os.path.exists(os.path.join(stat_root, sample)):
        os.makedirs(os.path.join(stat_root, sample))
        stat_pre = os.path.join(stat_root, sample, sample)
        c1 = "iteres stat {genome_size} {subfam_size} {rmsk} {bam} -o {pre} -T".format(
            genome_size=CHROM_SIZE,
            subfam_size=SUBFAMILY_SIZE,
            rmsk=RMSK,
            bam=bam,
            pre=stat_pre)
        call_sys([c1])
    if not os.path.exists(os.path.join(filter_root, sample)):
        os.makedirs(os.path.join(filter_root, sample))
        filter_pre = os.path.join(filter_root, sample, sample)
        #if unique mapping -N 0 means no difference to -N 2, however to test multiple mapping,it's necessary
        c2 = "iteres filter {genome_size} {subfam_size} {rmsk} {bam} -o {pre} -T".format(
            genome_size=CHROM_SIZE,
            subfam_size=SUBFAMILY_SIZE,
            rmsk=RMSK,
            bam=bam,
            pre=filter_pre)
        call_sys([c2])


def pipeline():
    mapping_root = "../2.Mapping/"
    stat_root = "./1.stat/"
    filter_root = "./2.filter/"
    bams = glob.glob(os.path.join(mapping_root, "*/*.bam"))
    bams.sort()
    Parallel(n_jobs=5)(delayed(iteres_stat_filter)(bam, stat_root, filter_root)
                       for bam in bams)


def getReps():
    reps = {}
    for line in open(REPEATS_LOCUS):
        line = line.split("\n")[0].split("\t")
        key = (line[0], line[1], line[2])
        reps[key] = "|".join(line)
    return reps


def getFPKMMatrix(reps, pre="TL_TEs"):
    fs = glob.glob("2.filter/*/*.loci")
    data = {}
    for f in fs:
        fn = f.split("/")[-2]
        print fn, " being processing."
        data[fn] = {}
        for line in open(f).read().split("\n")[1:]:
            line = line.split("\n")[0].split("\t")
            if len(line) < 3:
                continue
            key = (line[0], str(int(line[1]) + 1), line[2])
            if key not in reps:
                continue
            fpkm = float(line[8])
            if fpkm > 0:
                data[fn][reps[key]] = fpkm
    print "building matrix"
    data = pd.DataFrame(data)
    data = data.fillna(0)
    print data.shape
    """
    nis = [ t[0] for t in mat.itertuples() if np.sum(t[1:]) == 0]
    data = data.drop( nis )
    print data.shape
    """
    fn = pre + "_FPKM.txt"
    data.to_csv(fn, sep="\t", index_label="rep")
    print "finished!"


def getCountsMatrix(reps, pre="TL_TEs"):
    fs = glob.glob("2.filter/*/*.loci")
    data = {}
    for f in fs:
        fn = f.split("/")[-2].split("_")[0]
        print fn, " being processing."
        data[fn] = {}
        for line in open(f).read().split("\n")[1:]:
            line = line.split("\n")[0].split("\t")
            if len(line) < 3:
                continue
            key = (line[0], str(int(line[1]) + 1), line[2])
            if key not in reps:
                continue
            counts = int(line[7])
            if counts > 0:
                data[fn][reps[key]] = counts
    print "building matrix"
    data = pd.DataFrame(data)
    data = data.fillna(0)
    print data.shape
    fn = pre + "_counts.txt"
    data.to_csv(fn, sep="\t", index_label="rep")
    print "finished!"


def filterCountsMatrix(f,readsCut=10,sampleCut=10):
    mat = pd.read_table(f,sep="\t",index_col=0,header=0)
    for t in mat.itertuples():
        s = np.array(t[1:])
        s = s[s>readsCut]
        if len(s) < sampleCut:
            ns.append(t[0])
    mat = mat.drop(ns,axis=0)
    mat.to_csv(f.replace(".txt","_filter.txt"),sep="\t")
 

def main():
    #pipeline()
    #reps = getReps()
    #getCountsMatrix(reps)
    filterCountsMatrix("TL_TEs_counts.txt")


if __name__ == '__main__':
    start_time = datetime.now()
    main()
    elapsed = datetime.now() - start_time
    print "The process is done"
    print "Time used:", elapsed
