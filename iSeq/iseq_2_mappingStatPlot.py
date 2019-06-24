

#3rd
import brewer2mpl
import numpy as np
import pandas as pd

#global settings
import matplotlib as mpl
mpl.use("pdf")
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["figure.figsize"] = (4 * 0.8, 2.75 * 0.8)
mpl.rcParams["font.size"] = 8.0
from mpl_toolkits.axes_grid.inset_locator import inset_axes
colors = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors
import pylab
import seaborn as sns
sns.set_style("white")



def plotMmHg():
    mmat = pd.read_table("../5.mapping_mouse/MappingStat.txt",index_col=0,sep="\t")
    hmat = pd.read_table("../4.mapping_human/MappingStat.txt",index_col=0,sep="\t")
    mr = mmat["MappingRatio(%s)"]
    hr = hmat["MappingRatio(%s)"]
    stat = {"TotalReads":hmat["TotalReads"],"humanMappingRatio(%)": hr,"mouseMappingRatio(%)":mr}
    stat = pd.DataFrame(stat)
    stat.to_csv("MappingSummary.txt",sep="\t")
    ns = hr.index.intersection(mr.index)
    hr = hr[ns]
    mr = mr[ns]
    ncs = [n.split("_")[-1] for n in ns]
    ncs = {n:i for i, n in enumerate(list(set(ncs)))}
    cs = np.array([ ncs[ n.split("_")[-1] ] for n in ns])
    print(cs)
    fig, ax = pylab.subplots()
    for label, c in ncs.items():
        ns = np.where(cs==c)[0]
        #ax.scatter(hr[ns], mr[ns], color=colors[c], label=label,s=2)
        ax.scatter(hr[ns], mr[ns], color=colors[c], label=label,s=2)
    ax.legend(fontsize=6)
    ax.set_xlabel("human mapping ratio(%)")
    ax.set_ylabel("mouse mapping ratio(%)")
    pylab.tight_layout()
    pylab.savefig("1_human_mouse_mappingRatio.pdf")


def plotTotalMapRatio():
    mat = pd.read_table("/home/caoy7/caoy7/Projects/4.indexing/1.20190607_KZ1822/6.mapping_summary/MappingSummary.txt",index_col=0,sep="\t")
    x = mat["TotalReads"]
    y = mat["humanMappingRatio(%)"] + mat["mouseMappingRatio(%)"]
    fig,ax = pylab.subplots()
    ns = mat.index
    ncs = [n.split("_")[-1] for n in ns]
    ncs = {n:i for i, n in enumerate(list(set(ncs)))}
    cs = np.array([ ncs[ n.split("_")[-1] ] for n in ns])
    for label, c in ncs.items():
        ns = np.where(cs==c)[0]
        ax.scatter(x[ns], y[ns], color=colors[c], label=label,s=2)
    ax.legend(fontsize=6)
    ax.set_xscale("log")
    ax.set_xlabel("total reads")
    ax.set_ylabel("mappable ratio")
    pylab.tight_layout()
    pylab.savefig("2_total_mapratio.pdf")

def selSamples():
    """
    Select tnp pfv combinations according to cutoffs.
    """
    mat = pd.read_table("/home/caoy7/caoy7/Projects/4.indexing/1.20190607_KZ1822/6.mapping_summary/MappingSummary.txt",index_col=0,sep="\t")
    print(mat.shape)
    x = mat["TotalReads"]
    y = mat["humanMappingRatio(%)"] + mat["mouseMappingRatio(%)"]
    x = x[x>=10000] #more than 10**4 raw reads
    y = y[x.index] #more than 60% total mapping ratio
    y = y[y>=60]
    mat = mat.loc[y.index,]
    print(mat.shape)
    ns = mat.index
    ncs = [n.split("_")[-1] for n in ns]
    ncs = {n:i for i, n in enumerate(list(set(ncs)))}
    cs = np.array([ ncs[ n.split("_")[-1] ] for n in ns])
    ss = []
    hcut, lowcut = 70,10
    for label, c in ncs.items():
        ns = np.where(cs==c)[0]
        if label in ["Tnp01","Tnp02"]: #mouse specific 
            a = mat.ix[ns,"mouseMappingRatio(%)"]
            a = a[a>=hcut]
            b = mat.loc[a.index,"humanMappingRatio(%)"]
            b = b[b<=lowcut]
            print(label,b.shape)
            ss.extend( list(b.index))
        if label in ["Tnp05","Tnp06"]: #human specific 
            a = mat.ix[ns,"humanMappingRatio(%)"]
            a = a[a>=hcut]
            b = mat.loc[a.index,"mouseMappingRatio(%)"]
            b = b[b<=lowcut]
            print(label,b.shape)
            ss.extend( list(b.index))
        if label in ["Tnp03","Tnp04"]: #mouse and human mixed
            a = mat.ix[ns,"humanMappingRatio(%)"]
            a = a[ a>lowcut]
            a = a[ a<hcut ]
            b = mat.loc[a.index,"mouseMappingRatio(%)"]
            b = b[ b>lowcut ]
            b = b[ b<hcut ]
            print(label,b.shape)
            ss.extend( list(b.index) )
    mat = mat.loc[ss,]
    mat.to_csv("MappingSummary_filter.txt",sep="\t")


def plotMmHg2():
    mat = pd.read_table("./MappingSummary_filter.txt",index_col=0,sep="\t")
    mr = mat["mouseMappingRatio(%)"]
    hr = mat["humanMappingRatio(%)"]
    ncs = [n.split("_")[-1] for n in mat.index]
    ncs = {n:i for i, n in enumerate(list(set(ncs)))}
    cs = np.array([ ncs[ n.split("_")[-1] ] for n in mat.index])
    print(cs)
    fig, ax = pylab.subplots()
    for label, c in ncs.items():
        ns = np.where(cs==c)[0]
        label = label + ":%s"%len(ns)
        ax.scatter(hr[ns], mr[ns], color=colors[c], label=label,s=2)
    ax.legend(fontsize=6)
    ax.set_xlabel("human mapping ratio(%)")
    ax.set_ylabel("mouse mapping ratio(%)")
    pylab.tight_layout()
    pylab.savefig("3_human_mouse_mappingRatio_filtered.pdf")

 

def main():
    #plotMmHg()
    #plotTotalMapRatio()
    #selSamples()
    plotMmHg2()


if __name__ == "__main__":
    main()
