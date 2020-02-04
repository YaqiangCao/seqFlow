#sys
import random
from glob import glob
from copy import deepcopy
from datetime import datetime

#3rd
import umap
import hdbscan
import brewer2mpl
import numpy as np
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed
from sklearn import manifold
from sklearn.decomposition import PCA
from sklearn.cross_decomposition import CCA
from scipy.stats import ttest_ind as ttest

#global settings
import warnings
warnings.filterwarnings("ignore")

#3rd plotting setting
import brewer2mpl
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.use("pdf")
import seaborn as sns
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["figure.figsize"] = (4, 2.75)
mpl.rcParams["savefig.transparent"] = True
mpl.rcParams["savefig.bbox"] = "tight"
mpl.rcParams["font.size"] = 8.0
mpl.rcParams["font.sans-serif"] = "Arial"
mpl.rcParams["savefig.format"] = "pdf"
import pylab
sns.set_style("white")
from matplotlib.colors import ListedColormap
sns.set_style("white")
colors = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors
colors.extend(brewer2mpl.get_map('Set3', 'qualitative', 8).mpl_colors)


def getUmapPro(mat, pre, fout, n=100, min_dist=0):
    #mat = pd.read_csv(fin,sep="\t",index_col=0)
    Y = umap.UMAP(
        n_neighbors=n,
        n_components=2,
        min_dist=min_dist,
        init="random",
        metric="euclidean",
        random_state=123,
        n_epochs=1000,
    ).fit_transform(mat.values.T)
    ds = {"umap1": Y[:, 0], "umap2": Y[:, 1]}
    ds = pd.DataFrame(ds, index=mat.columns)

    labels = hdbscan.HDBSCAN(min_samples=20,
                             min_cluster_size=100).fit_predict(Y)
    clusters = pd.Series(labels, index=mat.columns)
    clusters = clusters[clusters >= 0]
    cs = (labels >= 0)

    fig, axs = pylab.subplots(1, 2, figsize=(8, 2.75))
    axs[0].scatter(Y[:, 0], Y[:, 1], c=colors[0], s=1)
    axs[0].set_title("raw %s UMAP projection" % pre)

    for c in set(labels):
        if c == -1:
            continue
        ps = np.where(labels == c)[0]
        axs[1].scatter(Y[ps, 0],
                       Y[ps, 1],
                       s=2,
                       c=colors[c],
                       label="cluster%s:%s cells" % (c, len(ps)))
    axs[1].set_title("UMAP projection and clustering")
    axs[1].legend(bbox_to_anchor=(1.02, 0.7), markerscale=5)
    pylab.savefig(fout + "_umap.pdf")

    clusters.to_csv(fout + "_cs.txt", sep="\t", header=False)
    mat = mat[clusters.index]
    mat.to_csv(fout + "_data.txt", sep="\t")
    ds = ds.loc[mat.columns, ]
    ds.to_csv(fout + "_umap.txt", sep="\t")


def showUmapGeneExp(expMat, umapMat, gene, pre):
    g = [g for g in expMat.index
         if g.split("|")[-1].upper() == gene.upper()][0]
    exp = expMat.loc[g, ]
    umapMat = umapMat.loc[exp.index, ]
    fig, ax = pylab.subplots(figsize=(4, 2.75))
    ax.scatter(umapMat["umap1"],
               umapMat["umap2"],
               c=exp,
               s=2,
               linewidths=0,
               cmap=plt.cm.Reds)
    ax.set_title(gene)
    ax.set_xlabel("UMAP-1")
    ax.set_ylabel("UMAP-2")
    pylab.savefig(pre + ".pdf")


def get():
    ds = pd.read_csv("../6.bigMat/allNorm.txt", index_col=0, sep="\t")
    print(ds.shape)

    wt = [c for c in ds.columns if "WT" in c]
    wtMat = ds[wt]
    getUmapPro(wtMat, "WT", "./WT")
    print("WT finished")

    n67 = [c for c in ds.columns if "N67" in c]
    n67Mat = ds[n67]
    getUmapPro(n67Mat, "N67", "./N67")
    print("N67 finished")

    ym = [c for c in ds.columns if "YM" in c]
    ymMat = ds[ym]
    getUmapPro(ymMat, "YM", "./YM")
    print("YM finished")


get()
