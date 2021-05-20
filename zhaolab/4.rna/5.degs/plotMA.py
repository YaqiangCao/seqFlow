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
    wt = mat["value_2"]
    ko = mat["value_1"]
    sa = wt[wt > 1].index
    sb = ko[ko > 1].index
    s = sa.union(sb)
    wt = wt[s]
    ko = ko[s]
    wt, ko = np.log2(wt), np.log2(ko)
    s = mat["significant"]
    s = s[s == "yes"]
    return wt, ko, s


def maPlot(wt, ko, pre, sig, mcut=np.log2(2)):
    fig, ax = pylab.subplots(figsize=(3.2, 2.2))
    sig = sig[sig > 0].index
    print(len(sig))
    #m = wt - ko
    m = ko - wt
    a = 0.5 * (wt + ko)
    up = m[m > mcut].index.intersection(sig)
    #up = a[up]
    #up = up[up>acut].index
    down = m[m < -mcut].index.intersection(sig)
    #down = a[down]
    #down = down[down>acut].index
    print(len(up), len(down))
    #all dots
    ax.scatter(a, m, color="gray", s=1, alpha=0.5)
    #up dots
    ax.scatter(a[up], m[up], color=colors[2], s=1, alpha=0.5)
    #down dots
    ax.scatter(a[down], m[down], color=colors[3], s=1, alpha=0.5)
    pylab.tight_layout()
    ax.set_xlabel("A, 1/2( log(GATA3 KO)+log(WT) )")
    ax.set_ylabel("M, log(GATA3 KO)-log(WT)")
    ax.axhline(y=0, linewidth=2, linestyle="--", color=colors[0])
    ax.axhline(y=mcut, linewidth=1, linestyle="--", color=colors[1])
    ax.axhline(y=-mcut, linewidth=1, linestyle="--", color=colors[1])
    ax.text(6, 2, "%s up genes" % len(up), color=colors[2])
    ax.text(6, -2, "%s down genes" % len(down), color=colors[3])
    # annotate the unchanged genes, between fc -1 and 1
    unc_up = m[m >= 0].index.difference(up)
    unc_down = m[m < 0].index.difference(down)
    ax.text(-2, 0.5, "%s genes" % len(unc_up), color="gray", fontsize=6)
    ax.text(-2, -0.5, "%s genes" % len(unc_down), color="gray", fontsize=6)
    #general settings
    ax.set_xlim([-2, 15])
    ax.set_ylim([-10, 15])
    ax.set_title(pre + ",%s genes" % len(m))
    pylab.savefig(pre + "_ma.pdf")


def main():
    for f in glob("*/gene_exp.diff"):
        n = f.split("/")[-2]
        wt, ko, sig = preMat(f)
        maPlot(wt, ko, n, sig)


main()
