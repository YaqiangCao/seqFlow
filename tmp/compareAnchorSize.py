import os, random, gzip
from glob import glob
from collections import Counter


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
        da = int(line[5]) -int(line[4])
        db = int(line[2]) - int(line[1])
        ds.append(da)
        ds.append(db)
    #ds = np.log2(np.array(ds))
    ds = np.array(ds)
    ds = ds[ds>0]
    return ds


def getdata(fs,fout):
    ds = {}
    i = 0
    for f in fs:
        n = f.split("/")[-1].split(".")[0].split("_")
        cell, tool = n[0],n[-1]
        print(cell,tool)
        dd = getds(f)
        for d in dd:
            ds[i] = {"tool":tool,"cell":cell,"size":d}
            i = i+1
    ds = pd.DataFrame(ds).T
    ds.to_csv(fout+".txt2",sep="\t")


def plotdata(f="all.txt"):
    mat = pd.read_table(f,index_col=0,sep="\t")
    #mat["size"] = mat["size"]/1000.0
    mat["size"] = np.log2(mat["size"])
    fig,ax = pylab.subplots(figsize=(4*0.8,2.75*0.8))
    sns.boxplot(x="tool",y="size",hue="cell",data=mat)
    pylab.savefig("all.pdf")



#getdata(glob("../1.sets/*.pgl"),"all")
plotdata()

