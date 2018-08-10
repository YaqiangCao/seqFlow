#!/usr/bin/env python2.7
#--coding:utf-8--
"""
getTree.py
2016-01-16: caculate distance matrix and prepare it for mega
"""

__author__ = "CAO Yaqiang"
__date__ = "2016-01-15"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import glob, os, time, string
from datetime import datetime

#3rd library
import pandas as pd
import numpy as np
from sklearn.metrics import pairwise_distances as pdist

#plot setting


def getdist(f):
    fout = f.split("/")[-1].replace(".txt", ".meg")
    mat = pd.read_table(f, index_col=0)
    dm = pdist(mat.T, metric='euclidean', n_jobs=10)
    dm = pd.DataFrame(dm, index=mat.columns, columns=mat.columns)
    dm = dm.mask(np.triu(np.ones(dm.shape)).astype(np.bool))
    dm.to_csv(fout, sep="\t", index=False, header=False)
    s = open(fout).read()
    cs = ["#%s" % c for c in mat.columns]
    cs = "\n".join(cs)
    ss = "#mega\n!TITLE  Genetic distance data;\n!Format DataType=distance;\n!Description\n    CYQ try;\n"
    ss = ss + cs + s
    with open(fout, "w") as f2:
        f2.write(ss)


def main():
    getdist("H3K4me1_TSS_sel.txt")


main()
