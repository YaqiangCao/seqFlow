#!/usr/bin/env python2.7
#--coding:utf-8--
"""
2019-05-21
"""

__author__ = "CAO Yaqiang"
__date__ = ""
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import os, time
from glob import glob
from datetime import datetime

#3rd library
import numpy as np
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed

#this
from utils import getLogger, callSys, cFlush

#global settings
#logger
#date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
#logger = getLogger(fn=os.getcwd() + "/" + date.strip() + "_" + os.path.basename(__file__) + ".log")


class PET(object):
    #cA is the center of left read
    __slots__ = [
        "chromA", "chromB", "startA", "startB", "endA", "endB", "strandA",
        "strandB", "cA", "cB", "distance", "cis"
    ]

    def __init__(self, d):
        """
        d is line = line.split( "\n" )[ 0 ].split( "\t" ) from BEDPE file 
        """
        self.chromA = d[0]
        self.startA = int(d[1])
        self.endA = int(d[2])
        self.strandA = d[8]
        self.chromB = d[3]
        self.startB = int(d[4])
        self.endB = int(d[5])
        self.strandB = d[9]
        self.cis = False
        if self.chromA == self.chromB:
            self.cis = True
            self.cA = (self.startA + self.endA) / 2
            self.cB = (self.startB + self.endB) / 2
            self.distance = abs(self.cB - self.cA)


def evaBedpe(f):
    total = 0  #total PETs
    cis = 0  #cis PETs
    close = 0  #distance < 150
    mid = 0  #distance > 150 < 1000
    far = 0  #distance > 20k
    fr = 0  #one mapped to postive and another mapped to negative strand
    ff = 0
    rr = 0
    for i, line in enumerate(open(f)):
        if i % 100000 == 0:
            cFlush("%s PETs processed from %s" % (i, f))
        line = line.split("\n")[0].split("\t")
        if len(line) < 6:
            continue
        total += 1
        pet = PET(line)
        if pet.cis:
            cis += 1
            if pet.distance <= 150:
                close += 1
            if 150 < pet.distance < 1000:
                mid += 1
            if pet.distance > 200000:
                far += 1
            if pet.strandA == "+" and pet.strandB == "-":
                fr += 1
            if pet.strandA == "-" and pet.strandB == "+":
                fr += 1
            if pet.strandA == "+" and pet.strandB == "+":
                ff += 1
            if pet.strandA == "-" and pet.strandB == "-":
                rr += 1
    sample = f.split("/")[-1].split(".")[0]
    s = {
        "totalPETs": total,
        "cisPETs": cis,
        "distance<=150bp": close,
        "150<=distance<=1000bp": mid,
        "distance>20kb": far,
        "forward/reverse": fr,
        "forward/forward": ff,
        "reverse/reverse": rr,
    }
    return sample, s


def main():
    fs = glob("*.bedpe")
    ds = {}
    for f in fs:
        sample, stat = evaBedpe(f)
        ds[sample] = stat
    ds = pd.DataFrame(ds).T
    ds.to_csv("bedpe_quality_stat.txt", sep="\t", index_label="sample")


if __name__ == '__main__':
    start_time = datetime.now()
    main()
    elapsed = datetime.now() - start_time
    print "The process is done"
    print "Time used:", elapsed
