#--coding:utf-8--
"""
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



def maPlot2(trt,con,pre,mcut=1,acut=2):
    fig, ax = pylab.subplots(figsize=(4,2.75))
    m = trt - con 
    a = 0.5 *  ( trt + con )
    up = m[m>mcut].index
    up = a[up]
    up = up[up>acut].index
    down = m[m<-mcut].index
    down = a[down]
    down = down[down>acut].index
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


def maPlot(mcut=1,acut=2):
    mat = pd.read_csv("KI_Rx1KI_Exp_filter.txt",index_col=0,sep="\t")
    cs = mat.columns
    neg = [ "KI_1","KI_2","KI_3"  ]
    pos = [ "Rx1KI_1","Rx1KI_2","Rx1KI_3" ]
    neg = mat[neg].mean(axis=1)
    pos = mat[pos].mean(axis=1)
    fig, ax = pylab.subplots(figsize=(4*0.8,2.75*0.8))
    m = pos - neg 
    a = 0.5 *  ( pos + neg )
    #all dots
    ax.scatter( a,m,color="gray",s=1,alpha=0.5)
    up = m[m>mcut].index
    down = m[m<-mcut].index
    #up dots
    ax.scatter(a[up],m[up],color=colors[2],s=1,alpha=1)
    #down dots
    ax.scatter(a[down],m[down],color=colors[3],s=1,alpha=1)

    ax.axhline(y=0,linewidth=1,linestyle="--",color=colors[0],alpha=0.5)
    ax.axhline(y=mcut,linewidth=1,linestyle="--",color=colors[1])
    ax.axhline(y=-mcut,linewidth=1,linestyle="--",color=colors[1])

    ax.text(3,4,"%s up genes"%len(up),color=colors[2])
    ax.text(3,-4,"%s down genes"%len(down),color=colors[3])
    # annotate the unchanged genes, between fc -1 and 1 
    unc_up = m[m>=0].index.difference( up )
    unc_down = m[m<0].index.difference( down )
    ax.text(9,1, "%s genes"%len(unc_up),color="gray")
    ax.text(9,-1, "%s genes"%len(unc_down),color="gray")
           

    ax.set_xlabel("A, 1/2( log(KI)+log(Rx1KI) )")
    ax.set_ylabel("M, log(Rx1KI) - log(KI)")
    pylab.savefig("expMA.pdf")
 
def main():
    maPlot()

main()
