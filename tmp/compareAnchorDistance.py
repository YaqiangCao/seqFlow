#sys
from glob import glob
import os, random, gzip

#3rd library
import HTSeq
from scipy.stats import wilcoxon 
from scipy.stats import mannwhitneyu
import matplotlib as mpl
from joblib import Parallel, delayed
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
import brewer2mpl
colors =  brewer2mpl.get_map( 'Set3', 'qualitative',8 ).mpl_colors




def getds(f):
    ds = []
    for line in open(f):
        line = line.split( "\n" )[ 0 ].split( "\t" )  
        d = (int(line[4]) + int(line[5]))/2 - (int(line[1])+int(line[2]))/2
        ds.append(abs(d))
    ds = np.array(ds)
    ds = ds[ds>1]
    return ds



def plot(pattern):
    fig,ax = pylab.subplots(figsize=(2*0.8,2.75*0.8))
    fs = glob("../1.sets/*%s*sorted"%pattern)
    fs.sort()
    for f in fs:
        d = getds(f)
        n = f.split("/")[-1].split(".")[0].split("_")[-1]
        sns.kdeplot(np.log2(d),ax=ax,label=n)
    ax.legend()
    ax.set_xlabel("distances between anchors, log2(bp)")
    ax.set_ylabel("density")
    pylab.savefig("%s.pdf"%pattern)
    

def plot2():
    fs = glob("../1.sets/*sorted")
    ds = {}
    i = 0
    for f in fs:
        print f
        diss = np.log2(getds(f))
        n = f.split("/")[-1].split(".")[0].split("_")
        tech = n[0]
        tool = n[1]
        for d in diss:
            ds[i] = {"dis":d,"tool":tool}
            i = i+1
    ds = pd.DataFrame(ds).T
    ds.to_csv("all_dis.txt",sep="\t")
    ds = pd.read_table("all_dis.txt",index_col=0)
    #ds["dis"] = np.log2(ds["dis"])
    fig,ax = pylab.subplots(figsize=(2.5,2.75))
    sns.violinplot(x="tool",y="dis",data=ds,palette="Set3",split=True,orient="v",saturation=0.5,ax=ax)
    ax.set_ylim([10,30])
    pylab.savefig("dis.pdf")

"""
data = {}
for f in glob("../1.sets/*.pgl"):
    n = f.split("/")[-1].split(".pgl")[0]
    ds = getds(f)
    print n,"\t", np.log2(ds.mean())
    data[n] = {"mean":np.log2(ds.mean()),"max":np.log2(ds.max())}
data = pd.DataFrame(data).T
data.to_csv("dis.txt",sep="\t")
"""


plot("_")
#plot("tech2")
#plot2()
