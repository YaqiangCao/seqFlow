#--coding:utf-8--
"""
getMA.py
2019-05-17: MA-plot based difference selection. M = logR - logG, A = 1/2(logR+logG)
"""

__author__ = "CAO Yaqiang"
__date__ = ""
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import os
from glob import glob

#3rd library
from tqdm import tqdm
import numpy as np
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
import pandas as pd
import brewer2mpl
colors = brewer2mpl.get_map('Set1', 'qualitative', 9).mpl_colors
colors.extend(brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors)
del colors[5]
del colors[10]


def preMat(f):
    mat = pd.read_table(f, index_col=0, sep="\t")
    cs = mat.columns
    wt = mat["value_1"]
    ko = mat["value_2"]
    sa = wt[wt > 1].index
    sb = ko[ko > 1].index
    s = sa.union(sb)
    wt = wt[s]
    ko = ko[s]
    wt, ko = np.log2(wt), np.log2(ko)
    s = wt - ko
    sa = s[s > 1]
    sb = s[s < -1]
    s = sa.index.union(sb)
    return wt, ko, s, mat.loc[s,]


def volcanoPlot(f, pcut=1e-3, fccut=1):
    n = f.split("/")[-2]
    
    mat = pd.read_csv(f, index_col=0, sep="\t")
    wt = mat["value_1"]
    ko = mat["value_2"]
    sa = wt[wt > 1].index
    sb = ko[ko > 1].index
    s = sa.union(sb)
    mat = mat.loc[s, ]

    p = mat["p_value"]
    fc = mat["log2(fold_change)"]
    sig = p[p < pcut]
    sig = fc[sig.index]
    up = sig[sig >= fccut].index
    down = sig[sig <= -fccut].index

    sig = up.union(down)
    fig, ax = pylab.subplots(figsize=(2,2.2))
    #all dots
    ax.scatter(fc, -np.log10(p), color="gray", s=1, alpha=0.5)
    #up dots
    ax.scatter(fc[up], -np.log10(p[up]), color=colors[2], s=2, alpha=0.5)
    #down dots
    ax.scatter(fc[down], -np.log10(p[down]), color=colors[3], s=2, alpha=0.5)

    ax.axhline(y=-np.log10(pcut), linewidth=1, linestyle="--", color=colors[0])
    ax.axvline(x=fccut, linewidth=1, linestyle="--", color=colors[1])
    ax.axvline(x=-fccut, linewidth=1, linestyle="--", color=colors[1])
    ax.text(2, 4.5, "%s up genes" % len(up), color=colors[2])
    ax.text(-10, 4.5, "%s down genes" % len(down), color=colors[3])
    ax.set_xlim([-10, 10])
    ax.set_ylim([0, 5])
    ax.set_xlabel("log2(fold change), KD/control")
    ax.set_ylabel("-log10(p-value)")
    pylab.savefig(n + "_volcano.pdf")
    #pylab.savefig(n + "_volcano.png")

    with open(n + "_up.list", "w") as fo:
        fo.write("\n".join([mat.loc[t, "gene"] for t in up]))
    with open(n + "_down.list", "w") as fo:
        fo.write("\n".join([mat.loc[t, "gene"] for t in down]))
    up = fc[up].sort_values(inplace=False,ascending=False)
    up.index = [ t + "|" + mat.loc[t,"gene"] for t in up.index]
    up.to_csv(n+"_up.rnk",sep="\t",header=None)
    down = fc[down].sort_values(inplace=False,ascending=True)
    down.index = [ t + "|" + mat.loc[t,"gene"] for t in down.index]
    down.to_csv(n+"_down.rnk",sep="\t",header=None)


def volcanoPlot2(f,ax, pcut=1e-3, fccut=1):
    n = f.split("/")[-2]
    
    mat = pd.read_csv(f, index_col=0, sep="\t")
    wt = mat["value_1"]
    ko = mat["value_2"]
    sa = wt[wt > 1].index
    sb = ko[ko > 1].index
    s = sa.union(sb)
    mat = mat.loc[s, ]

    p = mat["p_value"]
    fc = mat["log2(fold_change)"]
    sig = p[p < pcut]
    sig = fc[sig.index]
    up = sig[sig >= fccut].index
    down = sig[sig <= -fccut].index

    sig = up.union(down)
    #all dots
    ax.scatter(fc, -np.log10(p), color="gray", s=1, alpha=0.5)
    #up dots
    ax.scatter(fc[up], -np.log10(p[up]), color=colors[2], s=2, alpha=0.5)
    #down dots
    ax.scatter(fc[down], -np.log10(p[down]), color=colors[3], s=2, alpha=0.5)

    ax.axhline(y=-np.log10(pcut), linewidth=1, linestyle="--", color=colors[0])
    ax.axvline(x=fccut, linewidth=1, linestyle="--", color=colors[1])
    ax.axvline(x=-fccut, linewidth=1, linestyle="--", color=colors[1])
    ax.text(2, 4.5, "%s up genes" % len(up), color=colors[2])
    ax.text(-10, 4.5, "%s down genes" % len(down), color=colors[3])
    ax.set_xlim([-10, 10])
    ax.set_ylim([0, 5])
    ax.set_title(n)
    return ax


def main():
    """
    for f in glob("*/gene_exp.diff"):
        print(f)
        volcanoPlot(f)
    """
    fig,axs = pylab.subplots(2,3,figsize=(8,5),sharex=True,sharey=True)
    axs = axs.reshape(-1)
    for i,c in enumerate(["CTCF","RAD21","CTCF_RAD21","HCFC1","ZNF143","ZNF143_HCFC1"]):
        f = c+"/gene_exp.diff" 
        ax = axs[i]
        volcanoPlot2(f, ax)
        if i == 0:
            ax.set_xlabel("log2(fold change), KD/control")
            ax.set_ylabel("-log10(p-value)")
    #pylab.savefig("all_vol.png")
    pylab.savefig("all_vol.pdf")


main()
