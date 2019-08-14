#!/usr/bin/env python3.6
#--coding:utf-8 --
"""
getPeakOverlaps.py
2019-08-07:
"""

import os, subprocess
from glob import glob

#3rd
import matplotlib as mpl
mpl.use("pdf")
import pylab
import brewer2mpl
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm

#global settings
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["figure.figsize"] = (4, 2.75)
mpl.rcParams["font.size"] = 10.0
sns.set_style("white")
colors = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors


def getN(f):
    c = open(f).read().split("\n")
    c = [t for t in c if t.strip() != ""]
    c = len(c)
    return c


def getJI(fi, fj):
    """
    Get the overlap stats
    """
    ci = getN(fi)
    cj = getN(fj)
    cmd = "intersectBed -a %s -b %s -u > ./tmp.bed" % (fi, fj)
    subprocess.call(cmd, shell=True)
    cij = getN("./tmp.bed")
    ji = cij / 1.0 / (ci + cj - cij)
    return ci, cj, cij, ji


def pre(key="H1_H3K4me3"):
    fas = glob("../1.peaks/*%s*.bed" % key)
    if len(fas) == 0:
        return None
    fb = glob("../0.ENCODE/*%s*.bed" % key)
    if len(fb) == 0:
        return None
    else:
        fb = fb[0]
    ds = {}
    for f in fas:
        pre = f.split("/")[-1].split(".")[0]
        d = "_".join(pre.split("_")[:2])
        tool = "_".join(pre.split("_")[2:])
        ds[tool] = {"f": f, "ref": fb}
    return ds


def getStat(ds, pre):
    data = {}
    for tool in ds.keys():
        ci, cj, cij, ji = getJI(ds[tool]["f"], ds[tool]["ref"])
        data[tool] = {
            "Peaks": ci,
            "ENCODE_Peaks": cj,
            "OverlapPeaks": cij,
            "OverlapRatio": cij / 1.0 / ci,
            "JaccardIndex": ji,
        }
    data = pd.DataFrame(data).T
    data.to_csv(pre + ".txt", sep="\t")


def main():
    for cell in ["H1", "K562"]:
        for t in [
                "H3K4me1", "H3K4me2", "H3K4me3", "H3K27me3", "H3K27ac", "H2AZ"
        ]:
            key = cell + "_" + t
            print(key)
            ds = pre(key)
            if ds is not None:
                getStat(ds, key)


if __name__ == "__main__":
    main()
