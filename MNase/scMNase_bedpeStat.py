#!/usr/bin/env python2.7
#--coding:utf-8--
"""
scMNase_bedpeStat.py
2019-06-13: modified as add more stats, such as how many cn and sp PETs.
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


def getStat(f, dfilter=[80, 140, 180]):
    """
    Get the name, redundancy, total PETs, PETs distance mean, distance std for a bedpe file.
    """
    print(f)
    n = f.split("/")[-1].split(".bedpe")[0]
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
        if line[0] != line[3] or "_" in line[0]:
            continue
        t += 1
        s = min(int(line[1]), int(line[4]))
        e = max(int(line[2]), int(line[5]))
        r = (line[0], s, e)
        if r not in uniques:
            uniques.add(r)
            d = e - s
            ds.append(d)
            if d <= dfilter[0]:  #subnucleosome-sized particles
                sp += 1
            elif dfilter[1] <= d <= dfilter[2]:  #canonical nucleosomes
                cn += 1
            else:
                other += 1
    ds = np.array(ds)
    if t > 0:
        redundancy = 1.0 - len(uniques) / 1.0 / t
    else:
        redundancy = 0.0
    return n, t, len(uniques), redundancy, cn, sp, other, ds.mean(), ds.std()


def main():
    fs = glob("*.bedpe.gz")
    fs.extend(glob("*.bedpe"))
    data = Parallel(n_jobs=40)(delayed(getStat)(f) for f in fs)
    ds = {}
    for d in data:
        ds[d[0]] = {
            "totalMappedPETs": d[1],
            "uniquePETs": d[2],
            "redundancy": d[3],
            "canonicalNucleosomePETs": d[4],
            "subnucleosomeSizeParticlesPETs": d[5],
            "otherPETs": d[6],
            "fragmentLengthMean": d[7],
            "fragmentLengthStd": d[8]
        }
    ds = pd.DataFrame(ds).T
    ds.to_csv("stat.txt", sep="\t")


if __name__ == '__main__':
    start_time = datetime.now()
    main()
    elapsed = datetime.now() - start_time
    print "The process is done"
    print "Time used:", elapsed
