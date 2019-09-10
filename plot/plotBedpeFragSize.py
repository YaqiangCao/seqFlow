#!/usr/bin/env python2.7
#--coding:utf-8--
"""
bedpeFragSizePlot.py
"""

__author__ = "CAO Yaqiang"
__date__ = "2019-08-12"
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


def plotDis():
    fig,ax = pylab.subplots(figsize=(4*0.8,2.75*0.8))
    fs = glob("*.bedpe.gz")
    fs.sort()
    for f in fs:
        n = f.split(".bedpe.gz")[0]
        ds = []
        for i,line in enumerate(gzip.open(f)):
            if i % 10000 == 0:
                cFlush("%s read from %s"%(i,f))
            line = line.split("\n")[0].split("\t")
            s = min(int(line[1]), int(line[4]))
            e = max(int(line[2]), int(line[5]))
            d = e -s
            ds.append(d)
        ds = np.array(ds)
        sns.kdeplot(ds,label=n,ax=ax)
    #ax.axvline(x=80,linewidth=1,linestyle="--",color="gray")
    #ax.text(80,0.001,"80bp")
    #ax.axvline(x=120,linewidth=1,linestyle="--",color="gray")
    #ax.text(120,0.001,"120bp")
    #ax.axvline(x=140,linewidth=1,linestyle="--",color="gray")
    #ax.text(140,0.001,"140bp")
    #ax.axvline(x=180,linewidth=1,linestyle="--",color="gray")
    #ax.text(180,0.001,"180bp")
    #ax.set_xlim([0,300])
    ax.set_title("fragment length")
    ax.set_xlabel("length")
    ax.set_ylabel("density")
    pylab.tight_layout()
    pylab.savefig("bedpeFragSize.pdf")


def main():
    plotDis()


if __name__ == '__main__':
    start_time = datetime.now()
    main()
    elapsed = datetime.now() - start_time
    print "The process is done"
    print "Time used:", elapsed
