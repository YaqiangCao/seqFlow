#!/usr/bin/env python2.7
#--coding:utf-8 --
"""
"""

__author__ = "CAO Yaqiang"
__date__ = "2016-03-22"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import glob,os,copy
from collections import Counter

#3rd library
#plot setting
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
sns.set_style("whitegrid")
import brewer2mpl
colors = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors
colors.extend(brewer2mpl.get_map('Set3', 'qualitative', 8).mpl_colors)
#computating and stat setting.
import pandas as pd
import numpy as np
from scipy.stats import combine_pvalues as cps
from orangecontrib.bio.utils.stats import Hypergeometric as hy
from orangecontrib.bio.utils.stats import Binomial as bin
from orangecontrib.bio.utils.stats import FDR as fdr
from joblib import Parallel, delayed


def sgST(reps,repM,reptfs,tf):
    print "Estimating %s now"%tf
    b = bin()
    h = hy()
    N = repM.shape[0]
    n = reps.shape[0]
    reptfs = reptfs[ reptfs > 0 ].index
    m = reptfs.shape[0]
    reptfs = reptfs.intersection(reps)
    k = reptfs.shape[0]
    if m == 0 or k==0:
        print tf," no binding"
        return None
    bp = b.p_value(k, N, m, n)
    hp = h.p_value(k, N, m, n)
    cp = cps([hp, bp], method="stouffer")[1]
    es = float(k) / float(m) / float(n) * float(N)
    return m,k,hp,bp,cp,es


def repST(reps,repM,f,pre):
    mat = pd.read_table(f,index_col=0)
    mat = mat.loc[mat.index.intersection(repM.index)]
    ds = Parallel(n_jobs=20)(delayed(sgST)(reps,repM,copy.deepcopy(mat[c]),c) for c in mat.columns )
    data = {}
    for i,d in enumerate(ds):
        try:
            if type(d) == None or len(d) !=6:
                continue
        except:
            continue
        data[mat.columns[i]] = {
            "bg": d[0],
            "fg": d[1],
            "hy_p": d[2],
            "bin_p": d[3],
            "combinedP": d[4],
            "es": d[5]
            }
    data = pd.DataFrame(data).T
    ps = data["combinedP"].values
    qs = fdr(ps)
    data["FDR"] = qs
    data = data.sort("FDR")
    data.to_csv("%s_st.txt"%pre, sep="\t", index_label="TF")


def sgTest():
    bed = "../ESC_iPSC.bed"
    repF = "/picb/molsysbio/usr/caoyaqiang/1.Projects/17.Alu/1.GenomeReference/1.hg38/6.Repeats_Annotation/1.GenomicLocation/hg38Rep.txt"
    reps = pd.read_table(bed,index_col=3,header=None).index
    repM = pd.read_table(repF, index_col=0)
    #for t in ["h64","ecto","endo","mesendo","meso"]:
    #    f = "/picb/molsysbio/usr/caoyaqiang/1.Projects/17.Alu/7.Others/2.ENCODE/2.OtherPeaks/1.Nature_2015_ES_TF_Histone/2.hg38/3.Peaks2Repeats/%s.txt"%t
    #    repST(reps,repM,f,t)
    f = "/picb/molsysbio/usr/caoyaqiang/1.Projects/17.Alu/7.Others/4.Cistrome/1.hESCs/4.MPeaks2Repeats/hESC.txt"
    repST(reps,repM,f,"Cistrome_hESC")


def preMat(fs, qcut=1e-20,ecut=5):
    cs = ["h64","mesendo","endo","meso","ecto"]
    reps = set()
    mats = {}
    for f in fs:
        pre = f.split("_st")[0]
        mat = pd.read_table(f, index_col=0)
        mat.index = [ i.split("_")[1] for i in mat.index]
        mats[pre] = mat
        #qs = mat["FDR"]
        qs = mat["bin_p"]
        qs = qs[qs < qcut]
        es = mat.loc[qs.index,"es"]
        es = es[es>ecut]
        print qs.shape,es.shape,pre
        rs = set(es.index)
        reps.update(rs)
    reps = list(reps)
    mat_qs = {}
    mat_es = {}
    for pre, mat in mats.items():
        #mat_qs[pre] = mat.loc[reps, "FDR"]
        mat_qs[pre] = mat.loc[reps, "bin_p"]
        mat_es[pre] = mat.loc[reps, "es"]
    mat_qs = pd.DataFrame(mat_qs)[cs]
    mat_es = pd.DataFrame(mat_es)[cs]
    mat_qs[mat_qs < 1e-300] = 1e-300
    mat_qs[mat_qs > qcut] = 1
    mat_qs = 0.0 - np.log10(mat_qs)
    mat_qs = mat_qs.fillna(0.0)
    mat_es = mat_es.fillna(0.0)
    a = mat_es.sum(axis=1)
    a.sort(ascending=1)
    mat_qs = mat_qs.loc[a.index,:]
    #sort all samples
    mat_es = mat_es.loc[mat_qs.index,mat_qs.columns]
    mat_qs.to_csv("FDR.csv", sep="\t", index_label="rep")
    mat_es.to_csv("EnrichmentScore.csv", sep="\t", index_label="rep")


def bubleLegend(cmap="Blues"):
    sns.set_style("white")
    #plot legend, bubble size
    fig, ax = pylab.subplots()
    x = [0, 0, 0, 0, 0, 0]
    x2 = [1, 1, 1, 1, 1, 1]
    y = [1, 2, 3, 4, 5, 6]
    size = [10, 20, 40, 60, 80, 100]
    c = np.linspace(0.5, 3, 6)
    ax.scatter(x, y, s=size, marker="o", c=[0, 0, 0, 0, 0, 0])
    ax.scatter(x2, y, s=[40] * 6, marker="o", c=c, cmap=cmap)
    for i in xrange(len(x)):
        ax.text(x[i], y[i], ">%s" % size[i])
    for i in xrange(len(x2)):
        ax.text(x2[i], y[i], ">%s" % c[i])
    pylab.savefig("qs_bubbles_legend.pdf")


def enrichPlotBubble(cmap="Blues"):
    mat_qs = pd.read_table("FDR.csv", index_col=0).T
    #mat_qs[ mat_qs>100 ] = 100
    mat_es = pd.read_table("EnrichmentScore.csv", index_col=0).T
    x, y = [], []
    for i in xrange(mat_qs.shape[0]):
        for j in xrange(mat_qs.shape[1]):
            x.append(i)
            y.append(j)
    x, y = np.array(x) + 0.5, np.array(y) + 0.5
    size = np.ravel(mat_qs.values)
    cs = np.ravel(mat_es.values)
    fig, ax = pylab.subplots(figsize=(1,4 ))
    #fig, ax = pylab.subplots()
    ax.scatter(x, y, s=size/5, c=cs, cmap=cmap, marker="o")
    #ax.scatter(x, y, s=cs, c=size, cmap=cmap, marker="o")
    nx = np.arange(mat_qs.shape[0])
    ny = np.arange(mat_qs.shape[1])
    ax.set_xticks(nx)
    ax.set_yticks(ny)
    ax.set_xticklabels(map(str, mat_qs.index), rotation=90,size=10)
    ax.set_yticklabels(map(str, mat_qs.columns),size=10)
    ax.set_ylim([0, mat_qs.shape[1]])
    ax.set_xlim([0, mat_qs.shape[0]])
    pylab.savefig("qs_bubbles.pdf")


def ESplot(f="Cistrome_hESC_st.txt",cut=20,xcut=None,pre="ESC_TF"):
    sns.set_style("white")
    mat = pd.read_table(f,index_col=0)
    s = mat["es"]
    if xcut==None:
        xcut = s.mean()
    s = s[s>xcut]
    s.sort(ascending=0)
    #s.index = [i.split("_")[1].upper() for i in s.index]
    #s = s[0-cut:]
    s = s[:cut]
    s.sort(ascending=1)
    print s
    fig,ax = pylab.subplots(figsize=(0.6,2))
    x = np.arange(0,s.shape[0])
    ax.barh(x,s.values,0.6,color="gray")
    pylab.axvline(x=xcut,linewidth=1,color="red",linestyle="-")
    ax.set_yticks(x)
    ax.set_yticklabels(map(str,s.index),fontsize=8)
    ax.set_ylim([0,s.shape[0]])
    sns.despine()
    pylab.savefig("%s.pdf"%pre )


def main():
    #sgTest()
    #fs = ["h64_st.txt","meso_st.txt","mesendo_st.txt","endo_st.txt","ecto_st.txt"]
    #preMat(fs)
    #enrichPlotBubble()
    ESplot()
    #ESplot(f="h64_st.txt",pre="h64",xcut=5)


if __name__ == "__main__":
    main()
