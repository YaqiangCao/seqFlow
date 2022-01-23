#!/usr/bin/env python
#--coding:utf-8--
"""
"""

__author__ = "CAO Yaqiang"
__date__ = "2019-05-29"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

from glob import glob

#3rd library
import pylab
import pyBigWig
import numpy as np
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed

import matplotlib
matplotlib.use('pdf')
import brewer2mpl
colors = brewer2mpl.get_map('Set1', 'qualitative', 9).mpl_colors
colors.extend(brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors)
del colors[5]
del colors[10]





def getTss(bw,tssf,ext=2500,skipZeros=True):
    n = bw.split("/")[-1].split(".bw")[0]
    print(n)
    bw = pyBigWig.open(bw)
    x = range(-ext,ext+1)
    c = 0
    ys = np.zeros(len(x))
    for line in tqdm(open(tssf).read().split("\n")):
        line = line.split("\n")[0].split("\t")
        if len(line) < 3:
            continue
        chrom = line[0]
        center =int( (int(line[1]) + int(line[2] ))/2 ) 
        start = center -ext 
        end = center + ext+1
        if start < 0:
            continue
        try:
            s = np.array(bw.values(chrom, start, end))
        except:
            continue
        s = np.nan_to_num(s)
        if line[5] == "-":
            s = list(s)
            s.reverse()
            s = np.array(s)
        if skipZeros:
            if np.sum(s) > 0:
                ys = ys + s
                c += 1
        else:
            ys = ys + s
            c += 1
    ys = ys / c
    fig, ax = pylab.subplots()
    ax.plot(x, ys)
    ax.set_xlabel("Distance to TSS (bp)") 
    ax.set_ylabel("RPM")
    pylab.savefig(n+"_tss.pdf")


def main():
    tssf = "/home/caoy7/caoy7/Projects/0.Reference/1.hg38/2.annotations/gencode_v30_pcRNA_tss.bed"
    #tssf = "/home/caoy7/caoy7/Projects/0.Reference/2.mm10/2.annotations/gencode_vM21_pcRNA_tss.bed"
    fs = glob("../5.bdgBws/*.bw")
    Parallel(n_jobs=min(10,len(fs)))(delayed(getTss)(f,tssf) for f in fs)


if __name__ == '__main__':
    main()
