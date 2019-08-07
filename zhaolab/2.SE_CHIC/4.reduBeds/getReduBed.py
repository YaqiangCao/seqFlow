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


def thinBed(f):
    fout = f.split("/")[-1]
    if os.path.isfile(fout):
        return
    redus = set()
    print(f)
    with gzip.open(fout, "w") as fo:
        for i, line in enumerate(gzip.open(f)):
            line = line.split("\n")[0].split("\t")
            s = int(line[1])
            e = int(line[2])
            r = (line[0], s, e)
            if r in redus:
                continue
            else:
                redus.add(r)
            fo.write("\t".join(line) + "\n")


def main():
    fs = glob("../3.beds/*.bed.gz")
    data = Parallel(n_jobs=20)(delayed(thinBed)(f) for f in fs)


if __name__ == '__main__':
    start_time = datetime.now()
    main()
    elapsed = datetime.now() - start_time
    print "The process is done"
    print "Time used:", elapsed
