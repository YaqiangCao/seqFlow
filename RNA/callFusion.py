#!/usr/bin/env python2.7
#--coding:utf-8--
"""
callFusion.py
2018-07-16: based on STAR-Fusion v1.3.2
"""

__author__ = "CAO Yaqiang"
__date__ = "2018-07-05"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import os, time, commands, shutil
from glob import glob
from datetime import datetime

#3rd library
import matplotlib as mpl
mpl.use("pdf")
import seaborn as sns
import pylab
sns.set_style("white")
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["figure.figsize"] = (4 * 0.8, 2.75 * 0.8)
mpl.rcParams["figure.dpi"] = 100
mpl.rcParams["savefig.transparent"] = True
mpl.rcParams["savefig.bbox"] = "tight"
mpl.rcParams["font.size"] = 10.0
mpl.rcParams["font.sans-serif"] = "Arial"
mpl.rcParams["savefig.format"] = "pdf"
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.stats import ttest_rel

#my own
from utils import getlogger, call_sys

#global settings
#logger
date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
logger = getlogger(fn=os.getcwd() + "/" + date.strip() + "_" +
                   os.path.basename(__file__) + ".log")
#data

#starFusion="/picb/molsysbio2/caoyaqiang/2.WWX_RNA-seq/0.Tools/STAR-Fusion-STAR-Fusion-v1.4.0/STAR-Fusion"
starFusion = "/picb/molsysbio2/caoyaqiang/2.WWX_RNA-seq/0.Tools/STAR-Fusion-STAR-Fusion-v1.3.2/STAR-Fusion"
fusionlib = "/picb/molsysbio2/caoyaqiang/2.WWX_RNA-seq/1.Reference/4.STAR-Fusion/GRCh38_v27_CTAT_lib_Feb092018/ctat_genome_lib_build_dir/"
FA = "/picb/molsysbio/usr/caoyaqiang/1.Projects/18/1.GenomeReference/1.hg38/1.hg38_Sequence/hg38.fa"


def callFusion(jun):
    n = jun.split("/")[-2]
    cmd = "{starFusion} --genome_lib_dir {fusionlib} -J {jun} --output_dir {out}".format(
        starFusion=starFusion, fusionlib=fusionlib, jun=jun, out=n)
    call_sys([cmd])


def summaryFusion(cut=5, ncs=5):
    ds = {}
    for f in glob("*/star-fusion.fusion_predictions.abridged.tsv"):
        n = f.split("/")[-2]
        print(n)
        ds[n] = {}
        for i, line in enumerate(open(f)):
            line = line.split("\n")[0].split("\t")
            if i == 0:
                continue
            ds[n][line[0]] = int(line[1])
    #print(ds)
    ds = pd.DataFrame(ds)
    ds = ds.fillna(0)
    print(ds.shape)
    ns = []
    for t in ds.itertuples():
        nt = np.array(t[1:])
        nt = nt[nt >= cut]
        if len(nt) < ncs:
            ns.append(t[0])
    ds = ds.drop(ns)
    #ns = [c for c in ds.columns if np.sum(ds[c]) == 0]
    #ds = ds.drop(ns,axis=1)
    print(ds.shape)
    ds.to_csv("fusions.txt", sep="\t")


#based on ERCC
def normalizeFusion():
    ercc = pd.read_table("../3.ReadsCount/TL_ERCC_counts.txt.summary",
                         index_col=0).loc["Assigned", ]
    ercc.index = [i.split("/")[-1].split("_")[0] for i in ercc.index]
    mat = pd.read_table("fusions.txt", index_col=0)
    for c in mat.columns:
        mat[c] = mat[c] / ercc[c] * 10**6
    mat.to_csv("fusions_RPM.txt", sep="\t")


def plotFusion():
    mat = pd.read_table("fusions_RPM.txt", index_col=0)
    celllines = ["9810", "CCLP", "HT", "RBE"]
    lungControl = ["197111L", "197222L", "197333L"]
    lungCancer = ["197111T", "197222T", "197333T"]
    for t in mat.itertuples():
        cs = []
        iccs = []
        ds = {}
        for i, nt in enumerate(t):
            c = mat.columns[i - 1]
            if c in celllines:
                cat = "ICC_cellLine"
            else:
                if c in lungControl:
                    cat = "Lung_Control"
                else:
                    if c in lungCancer:
                        cat = "Lung_Cancer"
                    else:
                        if c.endswith("L"):
                            cat = "ICC_control"
                            cs.append(nt)
                        else:
                            cat = "ICC"
                            iccs.append(nt)
            ds[c] = {"RPM": nt, "cat": cat}
        ds = pd.DataFrame(ds).T
        fig, ax = pylab.subplots()
        p = ttest_rel(cs, iccs)[1]
        sns.stripplot(x="cat", y="RPM", data=ds, jitter=True)
        ax.set_ylim([0, 20])
        ax.set_title("%s\nICC vs control\npaired t-test p-value:%f" %
                     (t[0], p))
        labels = ax.get_xticklabels()
        ax.set_xticklabels(labels, rotation=45)
        print list(labels)
        pylab.savefig("%s.pdf" % t[0])


def main():
    #js = glob("../2.Mapping/*/*.junction")
    #Parallel( n_jobs=10 )( delayed( callFusion )( j ) for j in js )
    #summaryFusion()
    #normalizeFusion()
    plotFusion()


if __name__ == '__main__':
    start_time = datetime.now()
    main()
    elapsed = datetime.now() - start_time
    print "The process is done"
    print "Time used:", elapsed
