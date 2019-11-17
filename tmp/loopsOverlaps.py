#sys
import os, commands 
from glob import glob
#3rd
import numpy as np
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
mpl.rcParams[ "axes.titlesize" ] = "medium"
mpl.rcParams[ "axes.labelsize" ] = "medium"
import pylab
sns.set_style( "white" )
import brewer2mpl
colors = brewer2mpl.get_map('Set3', 'qualitative', 8).mpl_colors
import pandas as pd

pglt="/picb/molsysbio/usr/caoyaqiang/4.ENV/SG/sage-7.3-cyq/upstream/bio/pgltools/sh/pgltools"



def getStat():
    fs = glob("../../1.sets/GM12878*.pgl.sorted")
    fb = "./GM12878_CTCF_cLoops.pgl.sorted"
    data = {}
    for fa in fs:
        n = fa.split("/")[-1].split(".pgl.sorted")[0]
        dsa = len(open(fa).readlines())
        dsb = len(open(fb).readlines())
        cmd = "%s intersect -u -a %s -b %s | wc -l"%(pglt,fa,fb) 
        d = int(commands.getoutput(cmd))
        ji = float(d) / (dsa+dsb-d)
        ra = float(d) / dsa
        data[n] = {"La":dsa,"Lb":dsb,"Lab":d,"JI":ji,"ratio":ra}
    data = pd.DataFrame(data).T
    print data
    data.to_csv("stat.txt",sep="\t")

def getStatMango():
    fs = glob("../../1.sets/GM12878*.pgl.sorted")
    fb = "./GM12878_CTCF_mango.pgl.sorted"
    data = {}
    for fa in fs:
        n = fa.split("/")[-1].split(".pgl.sorted")[0]
        dsa = len(open(fa).readlines())
        dsb = len(open(fb).readlines())
        cmd = "%s intersect -u -a %s -b %s | wc -l"%(pglt,fa,fb) 
        d = int(commands.getoutput(cmd))
        ji = float(d) / (dsa+dsb-d)
        ra = float(d) / dsa
        data[n] = {"La":dsa,"Lb":dsb,"Lab":d,"JI":ji,"ratio":ra}
    data = pd.DataFrame(data).T
    print data
    data.to_csv("stat_usemango.txt",sep="\t")




def plotStat():
    data = pd.read_table("stat.txt",index_col=0)
    s = data["JI"]
    ds = {}
    i = 0
    for t in s.index:
        n = t.split("_")
        tool = n[-1]
        ds[i] = {"caller":'cLoops',"tool":tool,"Jaccard index":s[t]}
        i += 1
    data2 = pd.read_table("stat_usemango.txt",index_col=0)
    s = data2["JI"]
    for t in s.index:
        n = t.split("_")
        tool = n[-1]
        ds[i] = {"caller":'Mango',"tool":tool,"Jaccard index":s[t]}
        i += 1
    ds = pd.DataFrame(ds).T
    print ds
    fig,ax = pylab.subplots(figsize=(2.5,2.75))
    orders = ['Fit-Hi-C', 'GOTHiC', 'HOMER', 'HiCCUPS', 'cLoops', 'diffHic']
    ax = sns.barplot(x="tool", y="Jaccard index", hue='caller', data=ds,saturation=0.5, palette='Set3', order=orders)
    ax.set_title('GM12878\nHi-C overlap with CTCF ChIA-PET')
    ax.set_xlabel('')
    ax.set_ylabel("Jaccard index")
    ax.set_xticklabels(orders, rotation=35,ha="right")
    ax.set_xticks(range(len(orders)), minor=True)
    ax.tick_params(axis='both', direction='out', length=2, labelsize='medium', top=False, right=False)
    ax.legend( loc='best', mode=None, title='ChIA-PET loops\ncalled by', fontsize='small' )
    pylab.savefig("overlapwithCTCFChIAPET_together.pdf")

#getStat()
#getStatMango()
plotStat()
