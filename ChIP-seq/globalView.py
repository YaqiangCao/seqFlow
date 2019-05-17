#!/usr/bin/env python2.7
#--coding:utf-8--
"""
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
import numpy as np
import pandas as pd
import brewer2mpl
colors = brewer2mpl.get_map('Set1', 'qualitative', 9).mpl_colors
colors.extend(brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors)
del colors[5]
del colors[10]
from sklearn.decomposition import PCA, FastICA
from sklearn.cross_decomposition import PLSCanonical
from sklearn import manifold

#plotting settings


def pca_plot(mat, pre="test"):
    """
    Principle analysis and plot. 
    """
    pca = PCA(n_components=2)
    mat_r = pca.fit(mat.values.T).transform(mat.values.T)
    vs = pca.explained_variance_ratio_
    cs = {}
    i = 0
    f, ax = pylab.subplots()
    for i, label in enumerate(list(mat.columns)):
        #ax.scatter( mat_r[ i,0 ],mat_r[ i,1 ],label=label )
        ax.scatter(mat_r[i, 0], mat_r[i, 1], color="k")
        ax.text(mat_r[i, 0], mat_r[i, 1], label, fontsize=4)
    #for label in sorted( cs.keys(  ) ):
    #    ax.plot( 1,1,color=cs[label],label=label,markeredgecolor='none' )
    #leg = ax.legend( loc="upper left",fancybox=True,bbox_to_anchor=( 1,1 ) ,fontsize="x-small")
    #pylab.setp( leg.get_texts(  ) )
    ax.set_xlabel("PC1 (Explained Variance:%.3f)" % vs[0])
    ax.set_ylabel("PC2 (Explained Variance:%.3f)" % vs[1])
    #ax.set_xlim( [ -15,40 ] )
    pylab.savefig(pre + "_pca.pdf", bbox_inches='tight')


def mds_plot(mat, pre="test"):
    #mds = manifold.MDS(  )
    mds = manifold.TSNE()
    Y = mds.fit_transform(mat.values.T)
    cs = {}
    i = 0
    for label in list(mat.columns):
        nlabel = id2samples[ids[label]]
        if nlabel not in cs:
            cs[nlabel] = colors[i]
            i += 1
    f, ax = pylab.subplots()
    for i, label in enumerate(list(mat.columns)):
        ax.scatter(Y[i, 0], Y[i, 1], color=cs[id2samples[ids[label]]])
    for label in sorted(cs.keys()):
        ax.plot(1, 1, color=cs[label], label=label, markeredgecolor='none')
    leg = ax.legend(loc="upper left",
                    fancybox=True,
                    bbox_to_anchor=(1, 1),
                    fontsize="x-small")
    pylab.setp(leg.get_texts())
    ax.set_xlabel("MDS1")
    ax.set_ylabel("MDS2")
    #ax.set_xlim( [ -15,40 ] )
    pylab.savefig(pre + "_mds.pdf", dvi=1000, bbox_inches='tight')


def cor_heatmap(mat, pre="test"):
    fig, ax = pylab.subplots(figsize=(3, 3))
    cors = mat.corr()
    print(cors)
    mask = np.zeros_like(cors, dtype=np.bool)
    mask[np.triu_indices_from(mask, k=0)] = True
    sns.heatmap(cors,
                square=True,
                xticklabels=False,
                yticklabels=False,
                mask=mask,
                vmin=0.5)
    y = list(cors.index)
    y.reverse()
    ax.set_yticklabels(y, rotation=0, fontsize=5)
    pylab.savefig(pre + "_cor.pdf")


def main():
    #mat = pd.read_table("../10.getNoiseFilter/__TPM.txt",index_col=0)
    #pca_plot( mat,pre="all" )
    mat = pd.read_table("../10.getNoiseFilter/_-6_CD4_TPM.txt", index_col=0)
    cor_heatmap(mat, pre="_-6_CD4_TPM")
    #pca_plot( mat,pre="_-6_CD4_TPM.txt" )
    mat = pd.read_table("../10.getNoiseFilter/6_Th1_TPM.txt", index_col=0)
    cor_heatmap(mat, pre="_6_Th1_TPM")
    #pca_plot( mat,pre="6_Th1_TPM.txt" )
    #mat = pd.read_table("../10.getNoiseFilter/6_Th2_TPM.txt",index_col=0)
    #pca_plot( mat,pre="6_Th2_TPM.txt" )
    #mat = pd.read_table("../10.getNoiseFilter/24_Th1_TPM.txt",index_col=0)
    #pca_plot( mat,pre="24_Th1_TPM.txt" )
    #mat = pd.read_table("../10.getNoiseFilter/24_Th2_TPM.txt",index_col=0)
    #pca_plot( mat,pre="24_Th2_TPM.txt" )
    #cor_heatmap(mat,pre=pre)


main()
