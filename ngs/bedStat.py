#!/usr/bin/env python2.7
#--coding:utf-8--
"""
bedStat.py
Basic bedpe files quality control, as how many reads, redundancy, etc, work for paired-end ChIP-seq like data.
"""

__author__ = "CAO Yaqiang"
__date__ = "2019-05-23"
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
from joblib import Parallel, delayed

#this
#from utils import getLogger, callSys, PET
from utils import *

#global settings
#logger
date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
logger = getLogger(fn=os.getcwd() + "/" + date.strip() + "_" +
                   os.path.basename(__file__) + ".log")


def getStat(f):
    """
    Get the name, redundancy, total PETs, PETs distance mean, distance std for a bedpe file.
    """
    print(f)
    n = f.split("/")[-1].split(".bed")[0]
    if f.endswith(".gz"):
        of = gzip.open(f)
    else:
        of = open(f)
    #unique reads
    uniques = set()
    ds = []
    t = 0
    cn, sp, other = 0, 0, 0
    for i, line in enumerate(of):
        line = line.split("\n")[0].split("\t")
        if len(line) < 3:
            continue
        t += 1
        s = int(line[1])
        e = int(line[2])
        r = (line[0], s, e)
        if r not in uniques:
            uniques.add(r)
            ds.append(s-e)
    ds = np.array(ds)
    if t > 0:
        redundancy = 1.0 - len(uniques) / 1.0 / t
    else:
        redundancy = 0.0
    return n, t, len(uniques), redundancy, ds.mean(), ds.std()


def main(cpu=40):
    fs = glob("*.bed.gz")
    fs.extend(glob("*.bed"))
    cpu = min(cpu, len(fs))
    data = Parallel(n_jobs=cpu)(delayed(getStat)(f) for f in fs)
    ds = {}
    for d in data:
        ds[d[0]] = {
            "totalMappedReads": d[1],
            "uniqueReads": d[2],
            "redundancy": d[3],
        }
    ds = pd.DataFrame(ds).T
    ds.to_csv("bedStat.txt", sep="\t")


if __name__ == '__main__':
    start_time = datetime.now()
    main()
    elapsed = datetime.now() - start_time
    print "The process is done"
    print "Time used:", elapsed
