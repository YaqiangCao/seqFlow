#!/usr/bin/env python2.7
#--coding:utf-8 --
"""
map2others.py
Using bnMapper.py map the interesting repeats regions of human to other species.
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
    metaf = "/picb/molsysbio/usr/caoyaqiang/1.Projects/17.Alu/1.GenomeReference/1.hg38/19.LiftOverChains/order.txt"
    meta = pd.Series.from_csv(metaf, sep=":")
    chains = glob(
        "/picb/molsysbio/usr/caoyaqiang/1.Projects/17.Alu/1.GenomeReference/1.hg38/19.LiftOverChains/*.chain"
    )
    beds = glob("../../../14.Modules/3.JSDSelection/3.Modules/*.bed")
    cmds = []
    for chain in chains:
        fn = chain.split("/")[-1]
        if fn not in meta:
            continue
        s = meta[fn]
        try:
            os.mkdir(s)
        except:
            continue
        for bed in beds:
            fout = s + "/%s" % bed.split("/")[-1]
            cmd = "bnMapper.py -f BED4 -o {fout} --gap 20 --threshold 0.1 {fin} {chain}".format(
                fout=fout, fin=bed, chain=chain)
            #callSys([cmd])
            cmds.append(cmd)
    print cmds
    Parallel(n_jobs=14)(delayed(callSys)([c]) for c in cmds)


map2others()
