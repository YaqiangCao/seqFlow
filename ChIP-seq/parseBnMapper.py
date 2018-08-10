#!/usr/bin/env python2.7
#--coding:utf-8 --
"""
mappedParse.py
"""

__author__ = "CAO Yaqiang"
__date__ = ""
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#sys libray
import os
from glob import glob

#3rd library
import pandas as pd
from joblib import Parallel, delayed


def rs2bed(fn, rs):
    with open(fn, "w") as f:
        for r in rs:
            nr = r.split(":")[1].split("|")
            nnr = [nr[0], nr[1], nr[2], r]
            f.write("\t".join(nnr) + "\n")


def main():
    for s in glob("../2.match/*"):
        if not os.path.isdir(s):
            continue
        ns = s.split("/")[-1]
        try:
            os.mkdir(ns)
        except:
            pass
        fs = glob(s + "/*.bed")
        for f in fs:
            fout = "/".join(f.split("/")[-2:])
            if os.path.isfile(fout):
                continue
            rs = []
            for line in open(f):
                line = line.split("\n")[0].split("\t")
                hrep = line[3]
                orep = line[7]
                rs.append(hrep + ":" + orep)
            rs = list(set(rs))
            rs2bed(fout, rs)


if __name__ == "__main__":
    main()
