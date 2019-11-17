#!/usr/bin/env python2.7
#--coding:utf-8 --
"""
RepeatsEnrichment_bg.py
2015-08-31:Modified bg sets
2015-12-05:Modified as parse the sub-family, binomial test p-values added, combined p-values added.
2015-12-19: Modified as check if file already exists. 
"""

__author__ = "CAO Yaqiang"
__date__ = "2015-12-05"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import glob,os
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
from joblib import Parallel,delayed



def preFs(fes, fins):
    ds = {}
    for f in fes:
        fn = f.split("/")[-1].split("_")[0]
        ds[fn] = [f]
    for f in fins:
        fn = f.split("/")[-1].split("_")[0]
        if fn not in ds:
            continue
        ds[fn].append(f)
    return ds


def getRepFamily(f, reps):
    """
    Get family count for each input file.
    """
    mat = pd.read_table(f, index_col=0)
    rs = mat.index.intersection(reps.index)
    fams = list(reps[rs].values)
    cs = dict(Counter(fams))
    return cs


def repST(fe, fin, pre, reps):
    fout = pre + "_st.txt"
    if os.path.exists(fout):
        print "%s has been generated,return" % fout
        return
    fg = getRepFamily(fe, reps)
    bg = getRepFamily(fin, reps)
    N = sum(bg.values())
    n = sum(fg.values())
    data = {}
    for key in fg.keys():
        m = bg[key]
        k = fg[key]
        h = hy()
        b = bin()
        hp = h.p_value(k, N, m, n)
        bp = b.p_value(k, N, m, n)
        es = float(k) / float(m) / float(n) * float(N)
        cp = cps([hp, bp], method="stouffer")[1]
        data[key] = {
            "bg": m,
            "fg": k,
            "hy_p": hp,
            "bin_p": bp,
            "es": es,
            "combinedP": cp
        }
    data = pd.DataFrame(data).T
    ps = data["combinedP"].values
    qs = fdr(ps)
    data["FDR"] = qs
    data = data.sort("FDR")
    data.to_csv(fout, sep="\t", index_label="repFamily")


def filterRs(mat,cut=3):
    ns = []
    for t in mat.itertuples():
        s = np.array(t[1:])
        s = s[s>0]
        if len(s) < cut:
           ns.append(t[0]) 
    mat = mat.drop(ns)
    return mat

def sortCs(mat):
    f = "/picb/molsysbio/usr/caoyaqiang/1.Projects/17.Alu/3.NIH_RoadMap/1.ProcessedTagAlign/0.Meta/sel.txt"
    f2 = "/picb/molsysbio/usr/caoyaqiang/1.Projects/17.Alu/3.NIH_RoadMap/1.ProcessedTagAlign/0.Meta/sort.txt"
    ts = pd.read_table(f,index_col=0).index
    ids = {}
    for t in ts:
        t = t.split("|")
        ids[t[0]] = t[1]
    ids = pd.Series(ids)
    cs = []
    for line in open(f2):
        gs = line.split( "\n" )[ 0 ].split( ":" )[1].split(",")
        for g in gs:
            if g not in ids.index or ids[g] not in mat.columns:
                continue
            cs.append(g)
    cs = list(ids[cs].values)
    cs.reverse()
    mat = mat[cs]
    return mat
     

def preMat(fs=glob.glob("*.txt"), qcut=1e-10,cut=3):
    reps = set()
    mats = {}
    for f in fs:
        pre = f.split("_")[0]
        mat = pd.read_table(f, index_col=0)
        mats[pre] = mat
        qs = mat["FDR"]
        qs = qs[qs < qcut]
        rs = set(qs.index)
        reps.update(rs)
    reps = list(reps)
    mat_qs = {}
    mat_es = {}
    for pre, mat in mats.items():
        mat_qs[pre] = mat.loc[reps, "FDR"]
        mat_es[pre] = mat.loc[reps, "es"]
    mat_qs = pd.DataFrame(mat_qs)
    mat_es = pd.DataFrame(mat_es)
    mat_qs[mat_qs < 1e-100] = 1e-100
    mat_qs[mat_qs > qcut] = 1
    a = mat_qs.sum(axis=1)
    a.sort(ascending=True)
    mat_qs = 0.0 - np.log10(mat_qs)
    mat_qs = mat_qs.loc[a.index,:]
    #filter only one situation
    mat_qs = filterRs(mat_qs,cut=cut)
    #sort all samples
    mat_qs = sortCs(mat_qs)
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
    mat_qs = pd.read_table("FDR.csv", index_col=0)
    mat_es = pd.read_table("EnrichmentScore.csv", index_col=0)
    x, y = [], []
    for i in xrange(mat_qs.shape[0]):
        for j in xrange(mat_qs.shape[1]):
            x.append(i)
            y.append(j)
    x, y = np.array(x) + 0.5, np.array(y) + 0.5
    size = np.ravel(mat_qs.values)
    cs = np.ravel(mat_es.values)
    fig, ax = pylab.subplots(figsize=(4,6 ))
    ax.scatter(x, y, s=size/8, c=cs, cmap=cmap, marker="o")
    nx = np.arange(mat_qs.shape[0])
    ny = np.arange(mat_qs.shape[1])
    ax.set_xticks(nx)
    ax.set_yticks(ny)
    ax.set_xticklabels(map(str, mat_qs.index), rotation=90,size=8)
    ax.set_yticklabels(map(str, mat_qs.columns),size=5)
    ax.set_ylim([0, mat_qs.shape[1]])
    ax.set_xlim([0, mat_qs.shape[0]])
    pylab.savefig("qs_bubbles.pdf")


def main():
    """
    repf = "/picb/molsysbio/usr/caoyaqiang/1.Projects/17.Alu/1.GenomeReference/1.hg38/6.Repeats_Annotation/1.GenomicLocation/hg38Rep.txt"
    reps = pd.read_table(repf, index_col=0)
    reps = reps["family"]
    #enhancer/input sets
    fes = glob.glob("../../../6.ELR_PLR/1.Consolidated/5.activeEnhancers/*.txt")
    #inputs sets
    fins = glob.glob("../../../5.FPM/1.Consolidated/*.txt")
    ds = preFs(fes, fins)
    #get p-values
    #Parallel(n_jobs=len(ds))(delayed(repST)(value[0], value[1], key, reps)
    Parallel(n_jobs=10)(delayed(repST)(value[0], value[1], key, reps)
                             for key, value in ds.items())
    """
    preMat(  )
    enrichPlotBubble(  )
    bubleLegend(  )


if __name__ == "__main__":
    main()
