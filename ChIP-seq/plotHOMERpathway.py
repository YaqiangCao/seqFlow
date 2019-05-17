#!/usr/bin/env python3.6
#--coding:utf-8--
"""
plotHomerPathway.py
2016-02-27: modified plot type
2019-05-13: updated
"""

__author__ = "CAO Yaqiang"
__date__ = "2016-02-27"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import glob, os

#3rd library
import matplotlib as mpl
mpl.use("pdf")
import seaborn as sns
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["figure.figsize"] = (4, 3.5)
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
import brewer2mpl
colors = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors


def preprocess(fs=None,
               pcut=np.log10(1e-5),
               gcut=5,
               tcut=3000,
               root_dir="1.Parsed_GO"):
    if os.path.exists(root_dir) == False:
        os.mkdir(root_dir)
    for f in fs:
        nf = f.split("/")[-2:]
        sample, fname = nf[0], nf[1]
        sample = root_dir + "/" + sample
        fname = sample + "/" + fname
        if os.path.exists(sample) == False:
            os.mkdir(sample)
        mat = pd.read_table(f, index_col=0)
        s = mat["logP"]
        s = s[s < pcut]
        s = mat.loc[s.index, "Genes in Term"]
        s = s[s < tcut]
        s = mat.loc[s.index, "Target Genes in Term"]
        s = s[s > gcut]
        mat = mat.loc[s.index, :]
        if mat.shape[0] >= 1:
            mat.to_csv(fname, sep="\t")


def furthurparse(root_dir="1.Parsed_GO"):
    ds = glob.glob(root_dir + "/*")
    for d in ds:
        try:
            f = glob.glob(d + "/*.txt")[0]
        except:
            continue
        mat = pd.read_table(f, index_col=1)
        s = 0.0 - mat["logP"]
        s.to_csv(d + "/terms.terms", sep="\t")


def termBarh(numbercut=10,
             root_dir="1.Parsed_GO",
             plot_dir="2.Barh_Plots",
             x="-log10(p), HOMER Hypergeometric Test",
             title="Enriched GO Terms",
             width=0.6):
    if os.path.exists(plot_dir) == False:
        os.mkdir(plot_dir)
    fs = glob.glob(root_dir + "/*/*.terms")
    for f in fs:
        print(f)
        sample = f.split("/")[-2]
        s = pd.Series.from_csv(f, sep="\t")
        if s.shape[0] > numbercut:
            s = s[:numbercut]
        s = s[::-1]
        ns = s.shape[0] * 2
        xa = np.arange(0, ns, 2)
        xb = np.arange(1, ns, 2)
        fig, ax = pylab.subplots(figsize=(2.5, 4))
        ax.barh(xa, s.values, width, color=colors[1])
        for i, t in enumerate(s.index):
            ax.annotate(t, (0 + 0.1, xb[i]))
        pylab.xlabel(x)
        pylab.title(title)
        ax.set_yticks([])
        #ax.set_ylim([0, ns])
        sns.despine()
        pylab.savefig(plot_dir + "/" + sample + ".pdf")


root = "../../1.HOMER/go_*/"
fs = glob.glob(root + "biological_process.txt")
preprocess(fs=fs)
furthurparse()
termBarh()
