#!/usr/bin/env python2.7
#--coding:utf-8 --
"""
map2others.py
Using bnMapper.py map the interesting repeats regions of human to other species.
2016-05-04: specific map to mm10.
"""

__author__ = "CAO Yaqiang"
__date__ = "2016-01-08"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#sys libray
import os
from glob import glob

#3rd library
import pandas as pd
from joblib import Parallel, delayed


def callSys(cmds):
    """
    Call systematic commands without return.
    """
    for c in cmds:
        print c
        try:
            os.system(c)
        except:
            print "ERROR for %s" % c


def map2others():
    bedRoot = "/picb/molsysbio/usr/caoyaqiang/1.Projects/17.Alu/1.GenomeReference/1.hg38/24.OtherRMSKs/"
    cmds = []
    for s in glob("../1.map/*"):
        if not os.path.isdir(s):
            continue
        ns = s.split("/")[-1]
        try:
            os.mkdir(ns)
        except:
            continue
        fs = glob(s + "/*.bed")
        bed = bedRoot + ns + "/reps.bed"
        if not os.path.isfile(bed):
            print "ERROR for %s" % bed
        for f in fs:
            fout = "/".join(f.split("/")[-2:])
            if os.path.isfile(fout):
                continue
            #cmd = "intersectBed -a {fin} -b {ref} -wo > {fout}".format(fin=f,ref=bed,fout=fout)
            cmd = "/picb/molsysbio/usr/caoyaqiang/4.ENV/2.Bios/bedtools/bin.bk/intersectBed -a {fin} -b {ref} -wo > {fout}".format(
                fin=f, ref=bed, fout=fout)
            #cmd = "bedops --intersect {fin} {ref} > {fout}".format(fin=f,ref=bed,fout=fout)
            cmds.append(cmd)
    Parallel(n_jobs=30)(delayed(callSys)([c]) for c in cmds)


map2others()
