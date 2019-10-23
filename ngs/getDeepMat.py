#--coding:utf-8--
"""
Get the matrix from deeptools.
"""

__author__ = "CAO Yaqiang"
__date__ = ""
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import os,gzip
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



def getDeepMat(f,filterZero=True):
    ds = {}
    for i, line in enumerate(gzip.open(f)):
        if i == 0:
            continue
        line = line.split("\n")[0].split("\t")
        n = "|".join(line[:3])
        ss = list(map(float,line[6:]))
        if filterZero:
            if np.sum(ss) == 0.0:
                continue
        ss = pd.Series(ss,index=np.arange(len(ss)))
        ds[n] =ss
    ds = pd.DataFrame(ds).T
    ds.to_csv( f.replace(".txt.gz","_mat.txt"),sep="\t")


def main():
    getDeepMat("AID_TCF1_-6_CD4_Neg_Peaks_Norm_Profiles.txt.gz")

main()
