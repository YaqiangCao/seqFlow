#!/usr/bin/env python2.7
#--coding:utf-8--
"""
bed2SNPs.py
"""

import glob
from sklearn.neighbors import KDTree
from joblib import Parallel, delayed
#3rd library
import matplotlib as mpl
mpl.use("pdf")
import seaborn as sns
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["figure.figsize"] = (4, 2.75)
mpl.rcParams["figure.dpi"] = 100
mpl.rcParams["savefig.transparent"] = True
mpl.rcParams["savefig.bbox"] = "tight"
mpl.rcParams["font.size"] = 10.0
mpl.rcParams["font.sans-serif"] = "Arial"
mpl.rcParams["savefig.format"] = "pdf"
import pylab
sns.set_style("white")
import numpy as np
import pandas as pd
from scipy.stats import norm


def getloc(filename):
    inpf = open(filename)
    reploc = {}
    ss = {}
    for line in inpf:
        line = line.split("\n")[0].split("\t")
        #info = line[3].split('|')
        info = line
        chr = info[0]
        cent = (float(info[1]) + float(info[2])) / 2
        loc = [cent, 0]
        reploc.setdefault(chr, []).append(loc)
        ss.setdefault(chr, []).append(info[3])
    inpf.close()
    return reploc, ss


def getTree(filename):
    reploc, ss = getloc(filename)
    Trees = {}
    for chrom, loc in reploc.iteritems():
        Trees[chrom] = KDTree(loc)
    return Trees, reploc, ss


def getDist(queryfile, targetfile, fout):
    Trees, reploc, snps = getTree(targetfile)
    inpf = open(queryfile)
    pre = targetfile.split('/')[-1].split('_')[0]
    outf = open(fout, 'w')
    for line in inpf:
        info = line.split("\n")[0].split("\t")
        chr = info[0]
        if chr not in Trees:
            continue
        cent = (float(info[1]) + float(info[2])) / 2
        loc = [cent, 0]
        tree = Trees[chr]
        dis, ind = tree.query(loc, 1)
        dis = reploc[chr][ind[0][0]][0] - cent
        snp = snps[chr][ind[0][0]]
        line = [info[3], dis, snp]
        outf.writelines("\t".join(map(str, line)) + "\n")
    inpf.close()
    outf.close()


def plotDist(f, cut=1e4):
    fig, ax = pylab.subplots(figsize=(2, 2.75))
    d = pd.read_table(f, index_col=0, header=None)[1]
    n = d.shape[0]
    d = d[d < cut]
    d = d[d > -cut]
    ax = sns.distplot(d, bins=50, kde=False)
    ax.set_xlabel("Distance of ELR -> SNP,bp")
    ax.set_ylabel("Frequency")
    ax.set_xlim([-cut, cut])
    ax.set_title("%.2f ELRs" % (d.shape[0] / float(n) * 100))
    pylab.savefig("dist2SNPs.pdf")


if __name__ == '__main__':
    queryfile = "../1.beds/allELRs.bed"
    #targetfile = '/picb/molsysbio/usr/caoyaqiang/1.Projects/17.Alu/1.GenomeReference/1.hg38/9.SNPs_eQTL/4.GWAScatalogEBI/GWAScatalogEBI.bed'
    targetfile = '/picb/molsysbio/usr/caoyaqiang/1.Projects/17.Alu/1.GenomeReference/1.hg38/9.SNPs_eQTL/2.GRASP/2.Filtered/GRASP2eQTLs_hg38.bed'
    getDist(queryfile, targetfile, "allELRs_SNPs.txt")
    plotDist("allELRs_SNPs.txt")
