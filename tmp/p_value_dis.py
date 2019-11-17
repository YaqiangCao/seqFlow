#!/usr/bin/env python2.7
#--coding:utf-8--

"""
2015-07-23: PLS decomposition added, plot basic information modified.
"""


__author__="CAO Yaqiang"
__date__="2014-12-21"
__modified__="2015-07-23"
__email__="caoyaqiang0410@gmail.com"



#systematic library
import os
from glob import glob
from datetime import datetime

#3rd library
import matplotlib as mpl
mpl.use( "pdf" )
import seaborn as sns
mpl.rcParams[ "pdf.fonttype" ] = 42
mpl.rcParams[ "figure.figsize" ] = ( 4,2.75 )
mpl.rcParams[ "figure.dpi" ] = 100
mpl.rcParams[ "savefig.transparent" ] = True
mpl.rcParams[ "savefig.bbox" ] = "tight"
mpl.rcParams[ "font.size" ] = 10.0
mpl.rcParams[ "font.sans-serif" ] = "Arial"
mpl.rcParams[ "savefig.format" ] = "pdf"
import pylab
sns.set_style( "white" )
import numpy as np
import pandas as pd






def plot_PETs_dist(f):
    n = f.split("/")[-1].split(".")[0].split("_")
    cell = n[0]
    tool = n[-1]
    ds,ps = [],[]
    for line in open(f):
        line = line.split( "\n" )[ 0 ].split( "\t" ) 
        d = (int(line[4]) + int(line[5]))/2-(int(line[1])+int(line[2]))/2
        d = np.log2(np.abs(d))
        p = np.log2(int(float(line[-1])))
        ds.append(d)
        ps.append(p)
    pcc = np.corrcoef(ds,ps)[0,1]
    print cell,tool,pcc
    """
    df = pd.DataFrame({"distance":ds,"PETs":ps})
    fig,ax = pylab.subplots()
    df.plot(kind="hexbin",
        x="distance",
        y="PETs",
        ax=ax,
        colormap=pylab.cm.BuPu,
        reduce_C_function=np.max,
        bins="log")
    ax.set_title("PCC:%.3f"%pcc)
    ax.set_xlabel("distance between anchors, log(bp)")
    ax.set_ylabel("log2(rab)")
    pylab.savefig("%s_%s.pdf"%(tech,tool))
    """
    return cell,tool,pcc



i = 0
ds = {}
for f in glob("../1.sets/*.sorted"):
    cell,to,r = plot_PETs_dist(f)
    ds[i] = {"tool":to,"pcc":r,"cell":cell}
    i = i+1
ds = pd.DataFrame(ds).T
fig,ax = pylab.subplots(figsize=(2.5,2.75))
sns.barplot(x="tool",y="pcc",hue="cell",data=ds,saturation=0.5,palette="Set3")
pylab.savefig("all_PETs_dis.pdf")
