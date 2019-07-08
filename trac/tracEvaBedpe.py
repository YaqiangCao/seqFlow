#!/usr/bin/env python
#--coding:utf-8--
"""
2019-05-21
"""

__author__ = "CAO Yaqiang"
__date__ = ""
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import os, time,gzip
from glob import glob
from datetime import datetime

#3rd library
import numpy as np
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed

#this
from utils import getLogger, cFlush

#global settings
#logger
#date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
#logger = getLogger(fn=os.getcwd() + "/" + date.strip() + "_" + os.path.basename(__file__) + ".log")



class PET(object):
    """
    Paired-end tags / PETs object.
    """
    #cA is the center of left read
    __slots__ = [
        "chromA",
        "chromB",
        "startA",
        "startB",
        "endA",
        "endB",
        "strandA",
        "strandB",
        "cA",
        "cB",
        "distance",
        "cis",
        "length",
        "mid",
        "mapq",
        "start",
        "end",
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
        self.mapq = int(d[7])
        if self.chromA == self.chromB:
            self.cis = True
            #adjust the left end and right end to make sure left is alwasy small than right
            if self.startA + self.endA > self.startB + self.endB:
                self.startA, self.startB = self.startB, self.startA
                self.endA, self.endB = self.endB, self.endA
                self.strandA, self.strandB = self.strandB, self.strandA
            self.cA = (self.startA + self.endA) / 2
            self.cB = (self.startB + self.endB) / 2
            self.distance = abs(self.cB - self.cA)  #used for long-range PETs
            self.length = self.endB - self.startA  #fragment length, used for short-range PETs
            self.mid = (self.startA + self.endB) / 2  #middle of the fragment
            self.start = self.startA
            self.end = self.endB


def evaBedpe(f):
    """
    Qaulity control for the bedpe file.
    """
    total = 0  #total PETs
    cis = 0  #cis PETs
    close = 0  #distance < 150
    mid = 0  #distance > 150 < 1000
    far = 0  #distance > 1k < 10k
    extrame = 0 #distance > 20k
    fr = 0  #one mapped to postive and another mapped to negative strand
    ff = 0
    rr = 0
    uniques = 0
    trans = 0
    reds = set() #redundancy PETs
    chrMRatio = 0
    for i, line in enumerate(gzip.open(f)):
        if i % 100000 == 0:
            cFlush("%s PETs processed from %s" % (i, f))
        line = line.split("\n")[0].split("\t")
        if len(line) < 6:
            continue
        total += 1
        pet = PET(line)
        t = (pet.chromA, pet.chromB, pet.startA, pet.endA, pet.startB, pet.endB,pet.strandA,pet.strandB)
        if t not in reds:
            reds.add(t)
        else:
            continue
        if pet.cis !=True:
            trans += 1
            continue
        #only counting intra-chromosomal interaction PETs
        cis += 1
        if pet.distance <= 150:
            close += 1
        if 150 < pet.distance <= 1000:
            mid += 1
        if 1000 < pet.distance <= 10000:
            far +=1
        if pet.distance > 20000:
            extrame += 1
        if pet.strandA == "+" and pet.strandB == "-":
            fr += 1
        if pet.strandA == "-" and pet.strandB == "+":
            fr += 1
        if pet.strandA == "+" and pet.strandB == "+":
            ff += 1
        if pet.strandA == "-" and pet.strandB == "-":
            rr += 1
        if pet.chromA == "chrM":
            chrMRatio += 1

    sample = f.split("/")[-1].split(".")[0]
    s = {
        "highQPETs": total,
        "nonRedudantPETs": len(reds),
        "uniqueIntraPETs": cis,
        "uniqueInterPETs": trans,
        "distance<=150bp": close,
        "150<=distance<=1000bp": mid,
        "1kb<distance<=10k":far,
        "distance>20kb": extrame,
        "distance>20kb_ratio": extrame/1.0/cis,
        "forward/reverse": fr,
        "forward/forward": ff,
        "reverse/reverse": rr,
        "chrMRatio": chrMRatio/1.0/cis,
        "redundancy": 1- len(reds) / 1.0 / total
    }
    return sample, s


def main():
    fs = glob("./*.bedpe.gz")
    data = Parallel(n_jobs=len(fs))(delayed(evaBedpe)(f) for f in fs)
    ds = {}
    for t in data:
        sample, stat = t[0],t[1]
        ds[sample] = stat
    ds = pd.DataFrame(ds).T
    ds.to_csv("bedpe_quality_stat.txt", sep="\t", index_label="sample")


if __name__ == '__main__':
    start_time = datetime.now()
    main()
    elapsed = datetime.now() - start_time
    print("The process is done")
    print("Time used:", elapsed)
