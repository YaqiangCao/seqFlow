#!/usr/bin/env python2.7
#--coding:utf-8--
"""
getReduBedpe.py
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
    fout = f.split("/")[-1]
    if os.path.isfile(fout):
        return
    redus = set()
    print(f)
    with gzip.open(fout, "w") as fo:
        for i, line in enumerate(gzip.open(f)):
            line = line.split("\n")[0].split("\t")
            if line[0] != line[3] or "_" in line[0] or "_" in line[3]:
                continue
            #remove redudant PETs
            s = min(int(line[1]), int(line[4]))
            e = max(int(line[2]), int(line[5]))
            r = (line[0], s, e)
            if r in redus:
                continue
            else:
                redus.add(r)
            #shroten the name
            line[6] = str(i)
            fo.write("\t".join(line) + "\n")


def main():
    fs = glob("../3.bedpe/*.bedpe.gz")
    data = Parallel(n_jobs=20)(delayed(thinBedpe)(f) for f in fs)


if __name__ == '__main__':
    start_time = datetime.now()
    main()
    elapsed = datetime.now() - start_time
    print "The process is done"
    print "Time used:", elapsed
