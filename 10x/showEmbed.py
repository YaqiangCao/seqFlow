#sys
import random
from glob import glob
from collections import Counter
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
mpl.rcParams["savefig.transparent"] = True
mpl.rcParams["savefig.bbox"] = "tight"
mpl.rcParams["font.size"] = 8.0
mpl.rcParams["font.sans-serif"] = "Arial"
mpl.rcParams["savefig.format"] = "pdf"
import pylab
from matplotlib.colors import ListedColormap
sns.set_style("white")
colors = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors
colors.extend(brewer2mpl.get_map('Set3', 'qualitative', 8).mpl_colors)
colors.extend(brewer2mpl.get_map('Set1', 'qualitative', 8).mpl_colors)
del colors[6]
del colors[7]
del colors[7]


def plotEmbeding(mat, Y, title, xlabel, ylabel, pre):
    """
    """
    #prapare samples for different color
    cs = np.array(["_".join(c.split("_")[:-1]) for c in mat.columns])
    ncs = list(set(cs))
    ncs = {c: i for i, c in enumerate(ncs)}

    #plot
    f, ax = pylab.subplots(figsize=(3.2, 2.75 * 0.8))
    for c, i in ncs.items():
        ps = np.where(cs == c)[0]
        ax.scatter(Y[ps, 0], Y[ps, 1], color=colors[i], s=1, label=c)
    leg = ax.legend(
        loc="best",
        fancybox=True,
        markerscale=5,
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    pylab.savefig(pre + ".pdf")


def showPCA(mat, fout):
    pca = PCA(n_components=2)
    mat_r = pca.fit(mat.values.T).transform(mat.values.T)
    vs = pca.explained_variance_ratio_
    xlabel = "PC1 (Explained Variance:%.3f)" % vs[0]
    ylabel = "PC2 (Explained Variance:%.3f)" % vs[1]
    plotEmbeding(mat, mat_r, "PCA", xlabel, ylabel, fout + "_pca")
    mat_r = pd.DataFrame(mat_r, columns=["PC1", "PC2"], index=mat.columns)
    mat_r.to_csv(fout + "_pca.txt", sep="\t")


def showUmap(mat, fout):
    mat_r = umap.UMAP(n_neighbors=100,
                      n_components=2,
                      min_dist=0,
                      init="random",
                      metric="euclidean",
                      random_state=123,
                      n_epochs=1000).fit_transform(mat.values.T)
    xlabel = "UMAP-1"
    ylabel = "UMAP-2"
    plotEmbeding(mat, mat_r, "UMAP", xlabel, ylabel, fout + "_umap")
    mat_r = pd.DataFrame(mat_r, columns=["UMAP1", "UMAP2"], index=mat.columns)
    mat_r.to_csv(fout + "_umap.txt", sep="\t")


def showTsne(mat, fout):
    mat_r = manifold.TSNE(n_components=2,
                          perplexity=50,
                          n_iter=500,
                          init="pca").fit_transform(mat.values.T)
    xlabel = "tSNE-1"
    ylabel = "tSNE-2"
    plotEmbeding(mat, mat_r, "tSNE", xlabel, ylabel, fout + "_tsne")
    mat_r = pd.DataFrame(mat_r, columns=["tSNE1", "tSNE2"], index=mat.columns)
    mat_r.to_csv(fout + "_tsne.txt", sep="\t")


for f in glob("*_exp.txt"):
    print(f)
    mat = pd.read_csv(f, index_col=0, sep="\t")
    print(mat.shape)
    fout = f.split("/")[-1].split("_exp")[0]
    showPCA(mat, fout)
    showUmap(mat, fout)
    showTsne(mat, fout)
