#sys
import random
from glob import glob
from copy import deepcopy 
from datetime import datetime

#3rd
import numpy as np
import pandas as pd
#3rd plotting setting
import pylab
import brewer2mpl
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.use("pdf")
import seaborn as sns
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["savefig.transparent"] = True
mpl.rcParams["savefig.bbox"] = "tight"
mpl.rcParams["font.size"] = 8.0
mpl.rcParams["font.sans-serif"] = "Arial"
mpl.rcParams["savefig.format"] = "pdf"
import pylab
sns.set_style("white")
from matplotlib.colors import ListedColormap
sns.set_style("white")
colors = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors
colors.extend( brewer2mpl.get_map('Set3', 'qualitative', 8).mpl_colors )
colors.extend( brewer2mpl.get_map('Set1', 'qualitative', 8).mpl_colors )
del colors[6]
del colors[7]
del colors[7]


def readCs(f):
    ds = {}
    for line in open(f):
        line = line.split("\n")[0].split("\t")
        ds[ line[0] ] =  line[1] 
    ds = pd.Series(ds)
    return ds


def get():
    cs = readCs("../refinedCellTypes.txt")
    mat = pd.read_csv("../../12.cellCorr/all_scc.txt",sep="\t",index_col=0)
    print(mat)
    ds = {}
    for c in ["Spermatogonia","Spermatocyte","RoundSpermatid","Elongating","Sertoli"]:
        cells = cs[cs==c].index
        ncells = [ c for c in cells if c[0]=="N"]
        hcells = [ c for c in cells if c[0]=="H"]
        nmat = mat.loc[ncells,ncells].mean()
        hmat = mat.loc[hcells,hcells].mean()
        for i in nmat.index:
            ds[i] = {"diversity":1-nmat[i], "cellType":c, "condition":"Normal"}
        for i in hmat.index:
            ds[i] = {"diversity":1-hmat[i], "cellType":c, "condition":"Hypoxia"}
    ds = pd.DataFrame(ds).T
    ds.to_csv("diversity.txt",sep="\t")
        
 
def plot():
    ds = pd.read_csv("diversity.txt",sep="\t",index_col=0)
    ds["diversity"] = ds["diversity"].astype("float64")
    targetCells = ["Spermatogonia","Spermatocyte","RoundSpermatid","Elongating","Sertoli"]
    pat = sns.set_palette(sns.color_palette([colors[1],colors[2]]))
    fig, ax = pylab.subplots(figsize=(2,2.75*0.8))
    sns.violinplot(data=ds,x="diversity",y="cellType",hue="condition",orient="h",order=targetCells,split=True,inner="quartile")
    ax.set_ylabel("")
    pylab.savefig("diversity.pdf")

#get()
plot()
