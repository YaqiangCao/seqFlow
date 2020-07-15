import random
from glob import glob
from collections import Counter

#3rd
import numpy as np
np.seterr(all = 'ignore') 
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed
from scipy.stats import ttest_ind as ttest

#plot related
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
import pandas as pd
import brewer2mpl
colors = brewer2mpl.get_map('Set1', 'qualitative', 9).mpl_colors
colors.extend(brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors)
del colors[5]
del colors[10]


def sampleDiff(sa,sb,fcmax,fcmin,size=1000):
    """
    Sampling cell numbers to specific size and get the differential test results.
    @param sa: np.array
    @param sb: np.array
    @param size: int
    """
    #sampling 
    csa = list(range(len(sa)))
    csb = list(range(len(sb)))
    random.shuffle(csa)
    random.shuffle(csb)
    csa = csa[:size]
    csb = csb[:size]
    #estimating
    a = sa[csa]
    b = sb[csb]
    ea = a[a>0]
    eb = b[b>0]
    ea = float(len(ea)) / len(a)
    eb = float(len(eb)) / len(b)
    t,p = ttest(a,b)
    ma = np.mean(a)
    mb = np.mean(b)
    fc = np.log2(ma)-np.log2(mb)
    if fc > fcmax:
        fc = fcmax
    if fc < fcmin:
        fc = fcmin
    rs = {"fc":fc,"meanExp":ma,"expRatio":ea,"meanExpOthers":mb,"expRatioOthers":eb,"p-value":p,"expRatioDiff": ea-eb}
    rs = pd.Series(rs)
    return rs


def est(g,sa,sb,fcmax,fcmin,size=1000,repeats=100,fccut=1,pcut=1e-10,expr=0.1,fdrcut=0.01):
    """
    Estimating expression difference through sampling
    @param g: gene id
    @param sa: expression vector, pd.Series
    @param sb: expression vector, pd.Series
    @param fcmax: float, to replace np.inf
    @param fcmin: float, to replace -np.inf
    """
    sa = sa.values
    sb = sb.values
    data = {}
    for i in range(repeats):
        r = sampleDiff( sa,sb,fcmax,fcmin,size)
        data[i] = r 
    data = pd.DataFrame(data).T
    data = data.fillna(0.0)
    s = data.mean(axis=0)
    r = s.to_dict()
    if s["fc"] >= fccut and s["p-value"] <= pcut and s["expRatioDiff"] >= expr:
         sa = data["fc"]
         sa = sa[sa>=fccut]
         sb = data.loc[sa.index,"p-value"]
         sb = sb[sb<=pcut]
         sc = data.loc[sb.index,"expRatioDiff"]
         sc = sc[sc>=expr]
         fdr = 1 - float(len(sc))/data.shape[0]
         r["FDR"] = fdr
         if fdr < fdrcut:
             r["sig"] = 1
         else:
             r["sig"] = 0
    else:
         r["FDR"] = 1
         r["sig"] = 0
    return g, r
  


def checkDiff(mata,matb,fout):
    """
    @param mata: pd.DataFrame, matrix for cluster a
    @param matb: pd.DataFrame, matrix for other clusters
    @param fccut: float, log2 fold change cutoff
    @param pcut: float, t-test p-value cutoff
    @param expr: float, expressed cells 
    """
    print("Getting the fold change range")
    #get the fcmax and fcmin
    fcs = []
    for n in tqdm(mata.index):
        a = mata.loc[n,]
        b = matb.loc[n,]
        fc = np.log2(a.mean())-np.log2(b.mean())
        if -np.inf < fc < np.inf:
            fcs.append(fc)
    fcmax = np.max(fcs)
    fcmin = np.min(fcs)
    print("Carrying out permutation test")
    #ds = Parallel(n_jobs=60, backend="multiprocessing")(delayed(est)(sa[0],np.array(sa[1:]),np.array(sb[1:]),fcmax,fcmin) for sa,sb in tqdm(list(zip(mata.itertuples(),matb.itertuples()))))
    ds = Parallel(n_jobs=60, backend="multiprocessing")(delayed(est)(g,mata.loc[g,],matb.loc[g,],fcmax,fcmin) for g in tqdm(mata.index))
    print("Collecting permutation results")
    rs = {}
    for d in ds:
        rs[d[0]] = d[1]
    rs = pd.DataFrame(rs).T
    rs.to_csv(fout+".txt",sep="\t")
    sig = rs["sig"]
    sig = sig[sig>0]
    sig = rs.loc[sig.index,"fc"]
    sig = sig.sort_values(inplace=False,ascending=False)
    print(fout,"DEGs",len(sig))
    with open(fout+".list","w") as fo:
        gs = [g.split("|")[0] for g in sig.index]
        fo.write("\n".join(gs))
    return sig.index


def get(pre):
    cs = {}
    for line in open("../refinedCellTypes.txt"):
        line = line.split("\n")[0].split("\t")
        cs[ line[0] ] = line[1]
    cs = pd.Series(cs)
    ncs = pd.Series(Counter(cs.values))
    ncs = ncs.sort_values()
    ncs = list(ncs.index)
    ncs.sort()
    mat = pd.read_csv("../all_exp.txt",index_col=0,sep="\t")
    mat = mat[cs.index]
    gs = []
    cells = []
    for c in ncs:
        a = cs[cs==c].index
        print(c,len(a))
        b = cs.index.difference(a)
        rs = checkDiff( mat[a],mat[b],"cluster_%s"%(c))
        cells.extend(a)
        gs.extend(rs)
    mat = mat.loc[gs,cells]
    mat = mat.fillna(0)
    for t in mat.itertuples():
        nt = np.array(t[1:])
        nt = (nt - nt.mean())/nt.std()
        mat.loc[t[0],] = nt
    mat.to_csv("%s_clusters.txt"%pre,sep="\t")


get("all")
