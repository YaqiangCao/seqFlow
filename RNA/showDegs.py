import os
from glob import glob
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.stats import zscore
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
from matplotlib.colors import ListedColormap
import pandas as pd
import brewer2mpl
colors = brewer2mpl.get_map('Set1', 'qualitative', 9).mpl_colors
colors.extend(brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors)
del colors[5]
del colors[10]



def get(af,df, cell, time):
    mata = pd.read_csv(af,index_col=0,sep="\t")
    matb = pd.read_csv(df,index_col=0,sep="\t")
    mat = matb.loc[matb.index.intersection(mata.index),]
    
    wts = [ c for c in mat.columns if "WT_%s_%s"%(cell,time) in c]
    kos = [ c for c in mat.columns if "MLL4KO_%s_%s"%(cell,time) in c]
    cs = []
    cs.extend(wts)
    cs.extend(kos)
    
    mat = mat[cs]
    sa = mat[wts].mean(axis=1)
    sb = mat[kos].mean(axis=1)
    s = sa / sb 
    s = s.sort_values(inplace=False,ascending=False)
    mat = mat.loc[s.index,]
    mat.index = [ s.split("|")[1] for s in mat.index]
    cs = mat.columns
    cs = [c.replace("_%s_%s"%(cell,time),"") for c in cs]
    mat.columns = cs
    cs.sort()
    cs.reverse()
    mat = mat[cs]
    mat.to_csv("%s_%s_degs.txt"%(cell,time),sep="\t",index_label="gene")
    ds = {}
    for tmp in mat.itertuples():
        s = np.array(tmp[1:])
        s = zscore(s)
        ds[ tmp[0]] = s
    ds = pd.DataFrame(ds).T
    ds.columns = mat.columns
    ds = ds.loc[mat.index,]
    ds.to_csv("%s_%s_degs_zscore.txt"%(cell,time),sep="\t",index_label="gene")


def plot():
    fs = glob("*_zscore.txt")
    for f in fs:
        n = f.split("_degs_zscore.txt")[0]
        mat = pd.read_csv(f,index_col=0,sep="\t") 
        print(mat.columns)
        fig, ax = pylab.subplots(figsize=(1,2))
        cmap = sns.color_palette("RdBu_r", 11).as_hex()
        cmap[int(len(cmap) / 2)] = "#FFFFFF"
        cmap = ListedColormap(cmap)
        center = 0
        vmin = None
        ax = sns.heatmap( mat,
                          xticklabels=True,
                          yticklabels=False,
                          linewidths=0.0,
                          square=False,
                          center=0,
                          cmap = cmap,
                          ax=ax,
                          cbar_kws={"orientation": "horizontal"}
                          )
        ax.tick_params(axis='x', top=False, labeltop=True, labelbottom=False, direction='out')
        ax.set_xticklabels(mat.columns,fontsize=5,rotation=45,ha="left")
        """
        ax.tick_params(axis='y', top=False, labelright=True,labelleft=False, labelbottom=False, direction='out')
        ax.set_yticklabels(mat.index,fontsize=5,rotation=0)
        """
        cbar = ax.collections[0].colorbar
        cbar.ax.tick_params(labelsize=5)
        ax.set_ylabel("")
        wts = [ c for c in mat.columns if "WT" in c]
        kos = [ c for c in mat.columns if "KO" in c]
        sa = mat[wts].mean(axis=1)
        sb = mat[kos].mean(axis=1)
        s = sa - sb 
        ps = s[s>0]
        ns = s[s<0]
        ax.set_ylabel("%s up;%s down"%( len(ns),len(ps) ), fontsize=6)
        pylab.savefig("%s_degs.pdf"%n)
        with open("%s_up.list","w") as fo:
            fo.write("\n".join(ps.index))
        with open("%s_down.list","w") as fo:
            fo.write("\n".join(ns.index))


for cell in ["Th1","Th2"]:
    for time in ["24h","72h"]:
        fa = "../5.degs/%s_%s_degs.txt"%(cell,time)
        fb = "MLL4KO_filter.txt"
        get(fa,fb,cell,time)

plot()
