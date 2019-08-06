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
    mat = pd.read_table(f,index_col=0,sep="\t")
    cs = mat.columns
    #neg = [ c for c in cs if not ("AID" in c and "Pos" in c) ]
    neg = [ c for c in cs if ("AID" in c and "Neg" in c) ]
    pos = [ c for c in cs if ("AID" in c and "Pos" in c) ]
    neg = mat[neg].mean(axis=1)
    pos = mat[pos].mean(axis=1)
    #neg = np.log2( neg + 1.0)
    #pos = np.log2( pos + 1.0)
    neg = np.log2( neg )
    pos = np.log2( pos )
    return pos,neg


def maPlot(trt,con,pre,mcut=1,acut=2):
    fig, ax = pylab.subplots(figsize=(4,2.75))
    m = trt - con 
    a = 0.5 *  ( trt + con )
    up = m[m>mcut].index
    up = a[up]
    up = up[up>acut].index
    down = m[m<-mcut].index
    down = a[down]
    down = down[down>acut].index
    print(len(up),len(down))
    #ax.scatter( a,m,facecolors='white',edgecolors=colors[1],s=1,alpha=0.5 )
    #ax.scatter( a,m,color=colors[1],s=1,alpha=0.2 )
    #all dots
    ax.scatter( a,m,color="gray",s=1,alpha=0.5)
    #up dots
    ax.scatter(a[up],m[up],color=colors[2],s=1,alpha=0.5)
    #down dots
    ax.scatter(a[down],m[down],color=colors[3],s=1,alpha=0.5)
    pylab.tight_layout()
    ax.set_xlabel("A, 1/2( log(AID)+log(Control) )")
    ax.set_ylabel("M, log(AID) - log(Control)")
    ax.axhline(y=0,linewidth=2,linestyle="--",color=colors[0])
    ax.axhline(y=mcut,linewidth=1,linestyle="--",color=colors[1])
    ax.axhline(y=-mcut,linewidth=1,linestyle="--",color=colors[1])
    ax.text(6,2,"%s up peaks"%len(up),color=colors[2])
    ax.text(6,-2,"%s down peaks"%len(down),color=colors[3])
    # annotate the unchanged peaks, between fc -1 and 1
    unc_up = m[m>=0].index.difference( up )
    unc_down = m[m<0].index.difference( down )
    ax.text(0.2,0.5, "%s peaks"%len(unc_up),color="gray",fontsize=6)
    ax.text(0.2,-0.5, "%s peaks"%len(unc_down),color="gray",fontsize=6)
    #general settings 
    ax.set_xlim([0,10])
    ax.set_ylim([-4,4])
    ax.set_title(pre+",%s peaks"%len(m))
    pylab.savefig(pre+"_ma.pdf")


def main():
    for f in glob("../../10.getNoiseFilter/*.txt"):
        n = f.split("/")[-1].split("_TPM")[0]
        trt,con = preMat(f)
        maPlot( trt,con,n )

main()
