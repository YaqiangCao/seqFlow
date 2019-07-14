#!/usr/bin/env python2.7
#--coding:utf-8--
"""
globalCorrView.py
2019-05-31: updated t-SNE
2019-06-26: UMAP added
"""

__author__ = "CAO Yaqiang"
__date__ = "2016-01-14"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#sys
from glob import glob
from datetime import datetime

#3rd
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
import umap
import numpy as np
import pandas as pd
import brewer2mpl
colors = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors
#colors = brewer2mpl.get_map('Set1', 'qualitative', 9).mpl_colors
#colors.extend(brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors)
#del colors[5]
#del colors[10]
from sklearn.decomposition import PCA
from sklearn import manifold
from joblib import Parallel, delayed

#plotting settings


def plotEmbeding(mat, Y, title, xlabel, ylabel, pre):
    """
    Plot the embedding result, such from PCA, t-SNE. 
    mat: the matrix used
    Y: the component matrix, such as from PCA.
    """
    #prapare samples for different color
    cs = ["_".join(c.split("_")[:-1]) for c in mat.columns]
    cs = list(set(cs))
    cs = {c: i for i, c in enumerate(cs)}
    f, ax = pylab.subplots()
    for i, label in enumerate(list(mat.columns)):
        c = cs["_".join(label.split("_")[:-1])]
        ax.scatter(Y[i, 0], Y[i, 1], color=colors[c],s=5)
        #ax.text(Y[i,0],Y[i,1],label)
    for label, c in cs.items():
        ax.plot(1, 1, color=colors[c], label=label, markeredgecolor='none')
    #leg = ax.legend(loc="upper left",
    leg = ax.legend(loc="best",
                    fancybox=True,
                    #bbox_to_anchor=(1, 1),
                    fontsize="x-small")
    #pylab.setp(leg.get_texts())
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    pylab.savefig(pre + ".pdf")


def pca_plot(mat, pre="test"):
    """
    Principle analysis and plot. 
    """
    pca = PCA(n_components=2)
    mat_r = pca.fit(mat.values.T).transform(mat.values.T)
    vs = pca.explained_variance_ratio_
    xlabel = "PC1 (Explained Variance:%.3f)" % vs[0]
    ylabel = "PC2 (Explained Variance:%.3f)" % vs[1]
    plotEmbeding(mat, mat_r, "PCA", xlabel, ylabel, pre + "_pca")


def mds_plot(mat, pre="test"):
    mds = manifold.MDS(n_components=2)
    Y = mds.fit_transform(mat.values.T)
    plotEmbeding(mat, Y, "MDS", "MDS-1", "MDS-2", pre + "_mds")


def tsne_plot(mat, p=10, pre="test"):
    #tsne = manifold.TSNE(n_components=2, init="pca", perplexity=p,random_state=123)
    tsne = manifold.TSNE(n_components=2, perplexity=p)
    #embeding space normalization
    Y = tsne.fit_transform(mat.values.T)
    y_min,y_max = Y.min(0),Y.max(0)
    Y = (Y - y_min)/(y_max-y_min)
    xlabel = "t-SNE-1"
    ylabel = "t-SNE-2"
    #title = "perplexity=%s,init=PCA" % p
    title = "perplexity=%s" % p
    plotEmbeding(mat, Y, title, xlabel, ylabel, pre + "_tsne")


def umap_plot(mat,n=5,pre="test"):
    Y = umap.UMAP(n_neighbors=n,n_components=2,metric="euclidean",random_state=123,n_epochs=500).fit_transform(mat.values.T) 
    plotEmbeding(mat,Y,"UMAP","UMAP-1","UMAP-2",pre+"_umap")

def umap_plot_sup():
    """
    Supviszed UMAP dimension reduction.
    """
    pass


def main():
    ps = range(5,60,5) #parameters for tSNE
    ns = [5,10,15,20,30]
    for f in glob("*.txt"):
        print(f)
        mat = pd.read_table(f, index_col=0,sep="\t")
        n = f.split("/")[-1].split(".txt")[0]
        """
        pca_plot(mat, n)
        mds_plot(mat, n)
        Parallel(n_jobs=1)(delayed(tsne_plot)(mat,p,"%s_p_%s"%(n,p)) for p in ps)
        Parallel(n_jobs=1)(delayed(umap_plot)(mat,p,"%s_p_%s"%(n,p)) for p in ns)
        """
 


main()
