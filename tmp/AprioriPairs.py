#!/usr/bin/env python2.7
#--coding:utf-8 --
"""
AprioriTes.py
2016-08-24
2016-10-31: Modified the support and confidence definitation
"""

__author__ = "CAO Yaqiang"
__date__ = ""
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import os
from glob import glob

#3rd library
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
import pandas as pd
from joblib import Parallel,delayed
import brewer2mpl
colors =  brewer2mpl.get_map('Set1', 'qualitative', 9).mpl_colors    
colors.extend(brewer2mpl.get_map('Set3', 'qualitative', 8).mpl_colors)


def getReps(f,reps):
    rs = []
    terms = set()
    for line in open(f):
        line = line.split( "\n" )[ 0 ].split( "\t" ) 
        a,b = line[-2].strip(),line[-1].strip()
        rs.extend(a.split(","))
        rs.extend(b.split(","))
    reps = reps[rs]
    reps = reps.dropna()
    rs = []
    for i,line in enumerate(open(f)):
        line = line.split( "\n" )[ 0 ].split( "\t" ) 
        a,b = line[-2].strip(),line[-1].strip()
        a,b = reps[a.split(",")],reps[b.split(",")]
        a,b = a.dropna(),b.dropna()
        a,b = list(set(list(a.values))),list(set(list(b.values)))
        if len(a) == 0 and len(b)==0:
            continue
        rs.append([a,b])
        if len(a) == 0 or len(b)==0:
            continue
        for i in a:
            for j in b:
                #also count reverse pairs in 
                key = (i,j)
                terms.add(key)
                key = (j,i)
                terms.add(key)
    return rs,terms


def cusapir(f,reps):
    print f
    rs,terms = getReps(f,reps)
    ds = {}
    print f,len(rs),len(terms)
    for term in terms:
        ca,cb = 0,0
        a,b = term[0],term[1]
        for r in rs:
            if (a in r[0] and b in r[1]) or (a in r[1] and b in r[0]):
                ca += 1
            if (a in r[0]) or (a in r[1]):
                cb += 1
        ca,cb = float(ca)/len(rs),float(cb)/len(rs)
        #support
        s = ca
        #confidence
        c = ca/cb
        ds["->".join(term)] = {"co-occurance":ca,"occurance of #1":cb,"support":s,"confidence":c}
    pre = f.split("/")[-1].split(".")[0]
    ds = pd.DataFrame(ds).T
    ds.to_csv(pre+".txt",sep="\t")
    print "Estimation of %s finished."%pre


def plotsc(f,xcut=0.2,ycut=0.5):
    pre = f.split(".")[0]
    mat = pd.read_table(f,index_col=0)
    fig,ax = pylab.subplots()
    #support and confidence
    x,y = mat["support"],mat["confidence"]
    ax.scatter(x,y,c="k")
    for i in x.index:
        if x[i] > xcut and y[i] > ycut:
            ax.scatter(x[i],y[i],c="r")
            ax.text(x[i],y[i],i,fontsize=5)
    ax.set_xlabel("support")
    ax.set_ylabel("confidence")
    ax.axvline(x=xcut,linestyle="--",color="gray")
    ax.axhline(y=ycut,linestyle="--",color="gray")
    ax.set_xlim([0,1])
    ax.set_ylim([0,1])
    pylab.savefig(pre+".pdf")



def plotD(data,pre,scut=0.2,ccut=0.5):
    width = 0.6
    fig, ax = pylab.subplots(figsize=(2.5,1))
    rs = list(data.index)
    rs.reverse()
    data = data.loc[rs,]
    samples = list( data.index )
    ind = np.arange( len( data.index ) )
    s = data["support"].values
    ax.bar( ind,s,width,color="gray" )
    ax.set_ylabel( "support" )
    ax.set_xticks( ind + width )
    ax.set_xticklabels( samples,rotation=90,ha="right",fontsize=5  )
    ax.axhline(y=scut,linestyle="--",color="gray")
    ax2 = ax.twinx(  )
    c = data["confidence"].values
    ax2.plot( ind+width/2,c,color=colors[1] )
    for t in ax2.get_yticklabels(  ):
        t.set_color( colors[1] )
    ax2.axhline(y=ccut,linestyle="--",color=colors[1])
    ax2.set_ylabel( "confidence" )
    pylab.savefig( "Aprori_%s.pdf"%pre )




def plotSS(target="MIR->L2",xcut=0.2,ycut=0.5):
    ds = {} 
    for f in glob("*.txt"):
        pre = f.split(".")[0]
        mat = pd.read_table(f,index_col=0)
        if mat.loc[target,"support"] > xcut and mat.loc[target,"confidence"] > ycut:
            ds[pre] = {"support":mat.loc[target,"support"],"confidence":mat.loc[target,"confidence"]} 
    ds = pd.DataFrame(ds).T
    plotD(ds,target)



def main():
    """
    repf = "/picb/molsysbio/usr/caoyaqiang/1.Projects/17.Alu/1.GenomeReference/1.hg38/6.Repeats_Annotation/1.GenomicLocation/hg38Rep.txt"
    reps = pd.read_table(repf, index_col=0)["family"]
    #f = "../0.Sets/ENCODE_K562_H3K27ac.tsv"
    #cusapir(f,reps)
    fs = glob("../2.map2rep/*.pgl")
    Parallel( n_jobs=2 )( delayed( cusapir )( f,reps  ) for f in fs )
    map(plotsc,glob("*.txt"))
    """
    #plotSS("MIR->L2")
    plotSS("L2->MIR")


if __name__ == "__main__":
    main()
