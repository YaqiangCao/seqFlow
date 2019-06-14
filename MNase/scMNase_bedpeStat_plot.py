
import os, gzip
from glob import glob

#3rd
import HTSeq
import brewer2mpl
import numpy as np
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed

#global settings
import matplotlib as mpl
mpl.use("pdf")
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["figure.figsize"] = (4*0.8, 2.75*0.8)
mpl.rcParams["font.size"] = 10.0
from mpl_toolkits.axes_grid.inset_locator import inset_axes
colors = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors
import pylab
import seaborn as sns
sns.set_style("white")



def plotStat1():
    ds = pd.read_table("stat.txt",index_col=0,sep="\t")
    cs = ds.index
    ncs = ["_".join(c.split("_")[:-1]) for c in cs]
    nc = list(set(ncs))
    nc = { c:i for i,c in enumerate(nc) }
    ncs = np.array([nc[c] for c in ncs])

    #FRiP, total reads
    fig, ax = pylab.subplots()
    x = ds["uniquePETs"]
    y = ds["redundancy"]
    #for label, c in nc.items():
        #ns = np.where(ncs==c)[0]
        #ax.scatter(x[ns], y[ns], color=colors[c], label=label,s=2)
        #ax.scatter(1, 1, color=colors[c], label=label, s=0)
    ax.scatter(x, y, color=colors[0], label="scMNase_T",s=2)
    #ax.axhline(y=0.1,linewidth=1,linestyle="--",color="gray")
    #ax.axvline(x=10**4,linewidth=1,linestyle="--",color="gray")
    #ax.axvline(x=10**6,linewidth=1,linestyle="--",color="gray")
    ax.set_xscale("log")
    """
    inset_ax = inset_axes(ax, width="30%",height="30%",loc="upper left")
    bxs = []
    bys = []
    bts = []
    bcs = []
    i = 0
    ss = []
    for label, c in nc.items():
        ns = np.where(ncs==c)[0]
        ca = len(ns)
        bxs.append( i)
        i = i+0.2
        nx = x[ns]
        nx = nx[ nx > 10**4]
        nx = nx[nx<=10**6]
        ny = y[nx.index]
        ny = ny[ny>=0.1]
        cb = len(ny)
        ss.extend(ny.index)
        bys.extend( [ca,cb] )
        bcs.extend( [colors[c],colors[c]] )
        bts.append( label )
        bxs.append( i)
        i = i+0.6
    inset_ax.bar(bxs,bys,color=bcs,width=0.2)
    inset_ax.set_xticks([0,1,2])
    inset_ax.set_xticklabels(bts,fontsize=6,rotation=90)
    inset_ax.set_yticks([])
    for i, v in enumerate(bys):
        inset_ax.text(bxs[i], v+10, str(v),fontsize=6 )
    sns.despine(ax=inset_ax, top=True, right=True, left=True, bottom=False, offset=None, trim=False)
    #ax.legend(fontsize=8,markerscale=2)
    """
    ax.set_xlabel("uniuqe mapped reads (pairs)")
    ax.set_ylabel("redudancy")
    #ax.set_ylim([0,0.6])
    #pylab.tight_layout()
    pylab.savefig("1_readsRedudancy.pdf")
    #ds = ds.loc[ss,]
    #ds.to_csv("scMNase_stat_filter1.txt",sep="\t")
    

def plotStat2():
    ds = pd.read_table("stat.txt",index_col=0,sep="\t")
    cs = ds.index
    ncs = ["_".join(c.split("_")[:-1]) for c in cs]
    nc = list(set(ncs))
    nc = { c:i for i,c in enumerate(nc) }
    ncs = np.array([nc[c] for c in ncs])
    
    x = ds["canonicalNucleosomePETs"] / ds["uniquePETs"]
    y = ds["subnucleosomeSizeParticlesPETs"] / ds["uniquePETs"]
    fig, ax = pylab.subplots()
    #for label, c in nc.items():
    #    ns = np.where(ncs==c)[0]
    #    ax.scatter(x[ns], y[ns], color=colors[c], label=label,s=2)
    ax.scatter(x, y, color=colors[0], label="scMNase_T",s=2)
    ax.legend(fontsize=8,markerscale=2)
    ax.set_xlabel("Nucleosome reads ratio",fontsize=10)
    ax.set_ylabel("Subnucleosome particles reads ratio",fontsize=10)
    #pylab.tight_layout()
    pylab.savefig("2_cN_sP.pdf")
  

def plotStat3():
    ds = pd.read_table("stat.txt",index_col=0,sep="\t")
    cs = ds.index
    ncs = ["_".join(c.split("_")[:-1]) for c in cs]
    nc = list(set(ncs))
    nc = { c:i for i,c in enumerate(nc) }
    ncs = np.array([nc[c] for c in ncs])

    s = ds["fragmentLengthMean"]
    nds = []
    labels = []
    ccs = []
    fig, ax = pylab.subplots()
    """
    for label, c in nc.items():
        ns = np.where(ncs==c)[0]
        nds.append(s[ns])
        labels.append( label )
        ccs.append( colors[c])
    bx = ax.boxplot(nds,labels=labels,patch_artist=True)
    """
    bx = ax.boxplot(s)
    #for i, patch in enumerate(bx['boxes']):
    #    patch.set(facecolor=ccs[i])
    #sns.boxplot(data=nds,labels=labels)
    pylab.savefig("3_PETdistance.pdf")
   

def main():
    plotStat1()
    plotStat2()
    plotStat3()


if __name__ == "__main__":
    main()
