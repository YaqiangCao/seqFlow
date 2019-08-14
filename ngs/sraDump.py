#!/usr/bin/env python
#--coding:utf-8--

"""
sraDump.py
dump .sra into .fastq
"""


import os
from glob import glob
from joblib import Parallel, delayed

def dump( sra ):
    fn = sra.split("/")[-1].replace(".sra","")
    fq1 , fq2 = fn+"_1.fastq.gz", fn + "_2.fastq.gz"
    if os.path.isfile(fq1) and os.path.isfile(fq2):
        return
    cmd="fastq-dump --split-3 %s"%sra
    print(cmd)
    os.system( cmd )
    cmd = "rm %s"%sra
    print( cmd )
    os.system( cmd )

if __name__=="__main__":
    sras=glob.glob( "../1.SRA/*.sra" )
    Parallel(n_jobs=5)(delayed(dump)( sra ) for sra in sras)
