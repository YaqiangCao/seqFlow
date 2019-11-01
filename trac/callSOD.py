#!/usr/bin/env python
#--coding:utf-8 --
"""
findSOD.py
Find the self orgainized domain based on insulation score. 
Insulation score defination according to https://www.nature.com/articles/nature20158
X(x,s) = number of contacts between any pair of elements in the interval (x − s, x + s )
Modified as following
X(x,s) = (X(x,s) − X(x + s/2,s/2) − X(x − s/2, s/2))/(X(x+s/2,s/2)+X(x-s/2,s/2)
Directly comparing the reads than spanning a border and in the sperated up-stream and down-stream regions.
IS > 0.3 can be defined as SOD bordres, based on the observation on Trac-looping data.

2019-10-28: change the binsize and step to small number for Trac-looping data and increase the IS cutoffs for identificaiton of borders.
"""


__date__ = "2019-10-27"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#general library
import os, argparse
from glob import glob
from collections import Counter
from datetime import datetime
from argparse import RawTextHelpFormatter

#3rd library
import numpy as np
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed

#cLoops2
from cLoops2.ds import XY
from cLoops2.io import parseIxy



def findSodBorders(f,cut=0.08,close=500):
    """
    Find SOD.
    """
    print(f)
    ds = []
    for line in open(f):
        line = line.split("\n")[0].split("\t")
        if float( line[3] ) >= 1:
            continue
        ds.append( [ line[0],int(line[1]),int(line[2]),float(line[3])] )
 
    #find borders
    rs = []
    i = 0
    while i < len(ds):
        if ds[i][3] > cut:
            p = i
            s = ds[i][3]
            j = i+1
            while j < len(ds):
                if ds[j][3] > cut:
                    if ds[j][1] - ds[j-1][2] > close:
                        break
                    if ds[j][3] > s:
                        s = ds[j][3]
                        p = j
                    j = j+1
                else:
                    break
            #get the region
            rs.append([ds[i][0], ds[i][1], ds[j-1][2], s, ds[p][1],ds[p][2]])
            i = j 
        else:
            i = i + 1
    #stich borders
    nrs = []
    i = 0
    while i < len(rs):
        if i == len(rs) - 1:
            nrs.append(rs[i])
            break
        if rs[i + 1][1] - rs[i][2] > close:
            nrs.append(rs[i])
            i += 1
        else:
            p = i #record the sumit
            for j in range(i + 1, len(rs)):
                if rs[j][1] - rs[j - 1][2] > close:
                    break
                if rs[j][3] > rs[p][3]:
                    p = j
            j = j - 1
            r = [rs[i][0], rs[i][1], rs[j][2], rs[p][3],rs[p][4],rs[p][5] ]
            nrs.append(r)
            i = j + 1
    print(f,"finished",len(nrs))
    return nrs

def call(fs,fout,cut):
    print(fout,cut)
    ds = Parallel(n_jobs=30)(delayed(findSodBorders)( f, cut) for f in fs)
    rs = []
    for d in ds:
        rs.extend(d)
    with open(fout,"w") as fo:
        for r in rs:
            fo.write( "\t".join( map(str,r)) +"\n")

def main():
    for cut in np.arange(0.05,0.5,0.01):
    #for cut in [0.08,0.09]:
        fs = glob("./IS_bdgs/trac1_*.bdg")
        call(fs, "./beds/trac1_%.2f_border.bed"%cut,cut)
        fs = glob("./IS_bdgs/trac2_*.bdg")
        call(fs, "./beds/trac2_%.2f_border.bed"%cut,cut)

if __name__ == "__main__":
    main()
