#!/usr/bin/env python
#--coding:utf-8--
"""
"""

__author__ = "CAO Yaqiang"
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import os
import gzip

#3rd library
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator


def getRCSeq(seq):
    """
    Get the reverse complementary sequence
    """
    tab = str.maketrans("ACTG", "TGAC")
    return seq.translate(tab)[::-1]


def parseMeta(f):
    mat = pd.read_csv(f, index_col=0, sep="\t")
    ds = {}
    for t in mat.itertuples():
        ind = t[6] + "+" + getRCSeq(t[7])
        ind = ind.replace(" ", "")
        s = t[0] + "_" + t[11].replace(" ", "_").replace("-", "_")
        ds[ind] = s
    return ds


def get(fq1, fq2, metaf):
    ds = parseMeta(metaf)
    print(ds)
    fos = {}
    i = 0
    with gzip.open(fq1, "rt") as f1, gzip.open(fq2, "rt") as f2:
        for r1, r2 in zip(FastqGeneralIterator(f1), FastqGeneralIterator(f2)):
            i += 1
            if i % 100000 == 0:
                print("%s reads processed" % i)
            r1, r2 = list(r1), list(r2)
            rid = r1[0]
            irid = rid.split(":")[-1]
            if irid in ds:
                n = ds[irid]
                if irid not in fos:
                    fos[irid] = [
                        gzip.open(n + "_R1.fastq.gz", "wt"),
                        gzip.open(n + "_R2.fastq.gz", "wt")
                    ]
                fo1 = fos[irid][0]
                fo2 = fos[irid][1]
                fo1.write("@%s\n%s\n+\n%s\n" % (rid, r1[1], r1[2]))
                fo2.write("@%s\n%s\n+\n%s\n" % (rid, r2[1], r2[2]))


get("../0.fq/Undetermined_S0_L001_R1_001.fastq.gz",
    "../0.fq/Undetermined_S0_L001_R2_001.fastq.gz", "20221115_KZ2536.txt")
