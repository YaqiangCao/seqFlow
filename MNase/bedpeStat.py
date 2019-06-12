#!/usr/bin/env python2.7
#--coding:utf-8--
"""
bedpeStat.py
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


def thinBedpe(f):
    redus = set()
    with open(f + ".2", "w") as fo:
        for i, line in enumerate(open(f)):
            line = line.split("\n")[0].split("\t")
            #remove redudant PETs
            t = line[:6]
            t.extend([line[8], line[9]])
            t = tuple(t)
            if t in redus:
                continue
            redus.add(t)
            #remove the chr1_, chr2_dask
            if "_" in line[0] or "_" in line[3]:
                continue
            #shroten the name
            line[6] = str(i)
            fo.write("\t".join(line) + "\n")
    cmds = ["mv %s %s" % (f + ".2", f), "gzip %s" % f]
    callSys(cmds, logger)



def getStat(f):
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
    for i, line in enumerate(of):
        line = line.split("\n")[0].split("\t")
        if line[0] != line[3] or "_" in line[0]:
            continue
        t += 1
        s = min(int(line[1]), int(line[4]))
        e = max(int(line[2]), int(line[5]))
        r = (line[0],s,e)
        if r not in uniques:
            uniques.add(r)
            d = e -s
            ds.append(d)
    ds = np.array(ds)
    redundancy = 1.0 - len(uniques)/1.0/t
    return n, t, len(uniques), redundancy, ds.mean(), ds.std()



def main():
    fs = glob("*.bedpe.gz") 
    fs.extend(glob("*.bedpe"))
    data = Parallel(n_jobs=len(fs))(delayed(getStat)(f) for f in fs)
    ds = {}
    for d in data:
        ds[d[0]] = {"totalPETs": d[1],"uniquePETs":d[2],"redundancy":d[3], "distance mean": d[4], "distance std": d[5]}
    ds = pd.DataFrame(ds).T
    ds.to_csv("stat.txt", sep="\t")




if __name__ == '__main__':
    start_time = datetime.now()
    main()
    elapsed = datetime.now() - start_time
    print "The process is done"
    print "Time used:", elapsed
