from glob import glob
from collections import Counter
import numpy as np
import pandas as pd
import seaborn as sns
from cLoops2.settings import *


def get(pre):
    ds = {}
    for f in glob("../1.homer/*_ano.txt"):
        n = f.split("/")[-1].split("_ano.txt")[0]
        mat = pd.read_csv(f,index_col=0,sep="\t")
        ss = []
        for t in mat.itertuples():
            p = t[0]
            s = t[7]
            if "(" in s:
                s = s.split("(")[0].strip()
            ss.append( s )
        ss = pd.Series(Counter(ss))
        ds[n] = ss
    ds = pd.DataFrame(ds).T
    ds = ds.fillna( 0 )
    ds.to_csv(pre+"_stat_number.txt",sep="\t")
    nds = {}
    for t in ds.itertuples():
        nds[t[0]] = pd.Series(np.array(t[1:])/ np.sum( np.array(t[1:])),index=ds.columns)
    nds = pd.DataFrame(nds).T
    nds.to_csv(pre+"_stat_ratio.txt",sep="\t")
    ds = {}
    k = 0
    for i in nds.index:
        for j in nds.columns:
            ds[ k ] = {
                "location": j,
                "sample": i,
                "ratio": nds.loc[i,j]
                }
            k += 1
    ds = pd.DataFrame(ds).T
    ds.to_csv(pre+"_stat.txt",sep="\t")

def plot(pre):
    ds = pd.read_csv( pre +"_stat.txt",index_col=0,sep="\t")
    fig, ax = pylab.subplots(figsize=(3.2,2.2))
    sns.barplot( y="location",x="ratio", hue="sample", data=ds,orient="h")
    ax.set_ylabel( "" )
    ax.set_xlabel("peaks ratio")
    pylab.savefig(pre+"_stat.pdf")

get("ILC2vsTh2")
plot("ILC2vsTh2")
        
