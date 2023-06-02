#!/usr/bin/env python
#--coding:utf-8--
"""
"""

__author__ = "CAO Yaqiang"
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import os
import gzip
from glob import glob


#3rd library
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator



for f in glob("*R1*"):
    n = f.split("_R1")[0]
    print(n)
    ds = {}
    with gzip.open(f, "rt") as f:
        for i,r in enumerate(FastqGeneralIterator(f)):
            if i % 1000000 == 0:
                print("%s read from %s"%(i, n))
            r = list(r)
            rid = r[0]
            ind = rid.split(":")[-1]
            if ind not in ds:
                ds[ind] = 0
            ds[ind] += 1
    ds = pd.Series(ds)
    ds.to_csv(n+".txt",sep="\t",header=None)
