#!/usr/bin/env python2.7
#--coding:utf-8--
"""
globalCorrView.py
2019-05-31: updated t-SNE
2019-06-26: UMAP added
2019-07-18: updated plotEmbeding
"""

__author__ = "CAO Yaqiang"
__date__ = "2016-01-14"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#sys
import random
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
colors = brewer2mpl.get_map('Set1', 'qualitative', 9).mpl_colors
colors.extend(brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors)
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
    cs = np.array(["_".join(c.split("_")[:-1]) for c in mat.columns])
    ncs = list(set(cs))
    ncs = {c: i for i, c in enumerate(ncs)}
    f, ax = pylab.subplots()
    for c,i in ncs.items():
        ps = np.where( cs == c)[0]
        ax.scatter(Y[ps, 0], Y[ps, 1], color=colors[i],s=5,label=c)
    leg = ax.legend(
        loc="best",
        fancybox=True,
        #bbox_to_anchor=(1, 1),
        #fontsize="x-small"
    )
    #pylab.setp(leg.get_texts())
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    pylab.savefig(pre + ".pdf")


def plotClusterHeatmap(mat, pre):
    """
    Plot the hierarchal clustering result for samples (columns), do-not use this for maxtrix larger than 1000*1000, could be quite slow. 
    """
    cs = ["_".join(c.split("_")[:-1]) for c in mat.columns]
    cs = list(set(cs))
    cs = {c: i for i, c in enumerate(cs)}
    cs = [cs["_".join(c.split("_")[:-1])] for c in mat.columns]
    sample_colors = [colors[c] for c in cs]
    fig, ax = pylab.subplots(figsize=(4, 4))
    g=sns.clustermap(
        mat,
        #center=0.5,
        cmap="vlag",
        row_cluster=False,
        col_cluster=True,
        #z_score=1,
        #xticklabels=False,
        xticklabels=True,
        #yticklabels=False,
        yticklabels=True,
        #metric="manhattan",
        #row_colors=sample_colors,
        col_colors=sample_colors)
    #pylab.show()
    pylab.setp(g.ax_heatmap.get_xticklabels(), rotation=90)
    pylab.savefig(pre + "_cluster_heat.pdf")


def plotCorrHeatmap(mat, pre):
    """
    Plot the hierarchal clustering result for samples (columns), do-not use this for maxtrix larger than 1000*1000, could be quite slow. 
    """
    mat = mat.corr()
    cs = ["_".join(c.split("_")[:-1]) for c in mat.columns]
    cs = list(set(cs))
    cs = {c: i for i, c in enumerate(cs)}
    cs = [cs["_".join(c.split("_")[:-1])] for c in mat.columns]
    sample_colors = [colors[c] for c in cs]
    fig, ax = pylab.subplots(figsize=(4, 4))
    g = sns.clustermap(
        mat,
        #center=0.5,
        cmap="vlag",
        row_cluster=True,
        col_cluster=True,
        #z_score=1,
        xticklabels=True,
        yticklabels=True,
        row_colors=sample_colors,
        col_colors=sample_colors)
    pylab.setp(g.ax_heatmap.get_xticklabels(), rotation=90)
    #pylab.show()
    pylab.savefig(pre + "_corr_heat.pdf")


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


def rle_plot(mat,pre="test"):
    """
    Relative log expression. Mat should be raw counts matrix. 
    """
    mat = np.log2(mat+1)
    for t in mat.itertuples():
        s = np.array( t[1:] )
        m = np.median( s )
        s = s - m
        mat.loc[t[0],] = s
    fig, ax = pylab.subplots()
    sns.boxplot( data=mat,ax=ax,color=colors[0],fliersize=0 )
    ax.set_ylabel("Relative log expression")
    #pylab.xticks(rotation=90)
    ax.xaxis.set_tick_params(rotation=90)
    pylab.savefig(pre+"_rle.pdf")

def tsne_plot(mat, p=10, pre="test"):
    #tsne = manifold.TSNE(n_components=2, init="pca", perplexity=p,random_state=123)
    tsne = manifold.TSNE(n_components=2, perplexity=p)
    #embeding space normalization
    Y = tsne.fit_transform(mat.values.T)
    y_min, y_max = Y.min(0), Y.max(0)
    Y = (Y - y_min) / (y_max - y_min)
    xlabel = "t-SNE-1"
    ylabel = "t-SNE-2"
    #title = "perplexity=%s,init=PCA" % p
    title = "perplexity=%s" % p
    plotEmbeding(mat, Y, title, xlabel, ylabel, pre + "_tsne")


def umap_plot(mat, n=5, pre="test"):
    #Y = umap.UMAP(n_neighbors=n,n_components=2,metric="correlation",random_state=123,n_epochs=500).fit_transform(mat.values.T)
    Y = umap.UMAP(n_neighbors=n,
                  n_components=2,
                  metric="manhattan",
                  random_state=123,
                  n_epochs=500).fit_transform(mat.values.T)
    plotEmbeding(mat, Y, "UMAP", "UMAP-1", "UMAP-2", pre + "_umap")

def umap_shuffle_plot(mat, n=5, pre="test"):
    cs =list(mat.columns)
    random.shuffle( cs )
    mat = mat[cs]
    Y = umap.UMAP(n_neighbors=n,
                  n_components=2,
                  metric="manhattan",
                  random_state=123,
                  n_epochs=500).fit_transform(mat.values.T)
    plotEmbeding(mat, Y, "UMAP", "UMAP-1", "UMAP-2", pre + "_umap_shuffle")




def main():
    #ps = range(5, 60, 5)  #parameters for tSNE
    #ns = [5, 10, 15, 20, 30]
    for f in glob("*txt"):
        print(f)
        mat = pd.read_csv(f, index_col=0, sep="\t")
        n = f.split("/")[-1].split(".txt")[0]
        #rle_plot(mat,n)
        #pca_plot(mat, n)
        #mds_plot(mat, n)
        #Parallel(n_jobs=1)(delayed(tsne_plot)(mat,p,"%s_p_%s"%(n,p)) for p in ps)
        #Parallel(n_jobs=1)(delayed(umap_plot)(mat,p,"%s_p_%s"%(n,p)) for p in ns)
        #tsne_plot(mat,30,"%s_p_%s"%(n,30))
        #umap_plot(mat,30,"%s_p_%s"%(n,30))
        #umap_shuffle_plot(mat,30,"%s_p_%s"%(n,30))
        plotCorrHeatmap(mat, n)
        #plotClusterHeatmap(mat, n)


main()
