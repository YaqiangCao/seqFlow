#!/usr/bin/env python
#--coding:utf-8--
"""
2022-11-02: add enzyme option
2022-11-07: add support if a GSE has multiple SRR associated files
"""

#sys
import os
import time
import gzip
import subprocess
from glob import glob
from datetime import datetime

import click
import pandas as pd
from joblib import Parallel, delayed


def download(sample, srrs):
    """
    Download SRR files and save as .fastq.gz
    """
    if len(srrs) == 1:
        srr = srrs[0]
        c1 = "prefetch -p --max-size 200000000 %s" % srr
        c2 = "fastq-dump --split-3 ./%s" % srr
        c3 = "rm -fvr %s" % srr
        c4 = "mv %s_1.fastq 1.fastq/%s_R1.fastq" % (srr, sample)
        c5 = "mv %s_2.fastq 1.fastq/%s_R2.fastq" % (srr, sample)
        c6 = "ls 1.fastq/%s_R*.fastq | xargs -n 1 -P 2 gzip" % sample
        for c in [c1, c2, c3, c4, c5, c6]:
            print(c)
            os.system(c)
    else:
        os.mkdir( sample ) 
        for srr in srrs:
            c1 = "prefetch -p --max-size 200000000 %s -O %s" %(srr,sample)
            c2 = "fastq-dump --split-3 %s/%s -O %s"%(sample,srr,sample)
            c3 = "rm -fvr %s/%s"%(sample, srr)
            for c in [c1,c2,c3]:
                print(c)
                os.system(c)
        r1s = glob(sample+"/*_1.fastq")
        r2s = [ f.replace("_1.fastq","_2.fastq") for f in r1s]
        c1 = "cat %s > 1.fastq/%s_R1.fastq"%(" ".join(r1s), sample)
        c2 = "rm %s"%(" ".join(r1s))
        c3 = "cat %s > 1.fastq/%s_R2.fastq"%(" ".join(r2s), sample)
        c4 = "rm %s"%(" ".join(r2s))
        c5 = "rm -fvr %s"%sample
        c6 = "ls 1.fastq/%s_R*.fastq | xargs -n 1 -P 2 gzip" % sample
        for c in [c1, c2, c3, c4, c5, c6]:
            print(c)
            os.system(c)



def run(sample, srrs,):
    """
    run the process for 1 sample
    """
    print(sample, srrs)
    fq1 = "1.fastq/%s_R1.fastq.gz" % sample
    fq2 = "1.fastq/%s_R2.fastq.gz" % sample
    if not (os.path.isfile(fq1) and os.path.isfile(fq2)):
        download(sample, srrs)
   

@click.command()
@click.option(
    "-cpu",
    required=True,
    help="CPUs for parallel samples.",
    default=10,
    type=int,
)
def main(cpu):
    wd = os.path.dirname(os.path.realpath(__file__))
    print("%s/run.py" % (wd))
    logf = "run.log"

    for d in ["1.fastq"]:
        if not os.path.exists(d):
            os.mkdir(d)
    #parse the meta information
    mat = pd.read_table("SraRunInfo.csv", index_col=0, sep=",")
    metamat = pd.read_table("sra_result.csv", index_col=0, sep=",")
    ds = {}
    for t in metamat.itertuples():
        srx = t[0]
        #sample = t[2]
        #sample = t[1].split(":")[0].strip()
        sample = metamat.loc[srx,"Experiment Title"].split(";")[0].strip().replace(": ","_")
        ss = mat["Experiment"]
        ss = ss[ss == srx].index
        ds[sample] = list(ss)
    Parallel(n_jobs=cpu)(delayed(run)(sample, srrs) for sample, srrs in ds.items())
    

if __name__ == "__main__":
    startTime = datetime.now()
    main()
    elapsed = datetime.now() - startTime
    print("The process is done.")
    print("Time used:", elapsed)
