from glob import glob
from collections import Counter
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.stats import ttest_ind as ttest
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



def checkDiff(mata,matb,fout,fccut=1,pcut=1e-5,expr=0.2):
    """
    @param mata: pd.DataFrame, matrix for cluster a
    @param matb: pd.DataFrame, matrix for other clusters
    @param fccut: float, log2 fold change cutoff
    @param pcut: float, t-test p-value cutoff
    @param expr: float, expressed cells 
    """
    s = mata.sum(axis=1)
    s = s[s>=s.median()]
    ns = s.index
    rs = {}
    ts = []
    for n in tqdm(mata.index):
        a = mata.loc[n,]
        b = matb.loc[n,]
        ea = a[a>0]
        eb = b[b>0]
        ea = float(len(ea)) / len(a)
        eb = float(len(eb)) / len(b)
        t,p = ttest(a,b)
        ma = np.mean(a)
        mb = np.mean(b)
        fc = np.log2(ma)-np.log2(mb)
        rs[n] = {"fc":fc,"meanExp":ma,"expRatio":ea,"meanExpOthers":mb,"expRatioOthers":eb,"p-value":p,"expRatioDiff": ea-eb}
        if fc > fccut and p < pcut and ea >= expr and n in ns:
            rs[n]["sig"] = 1
            ts.append( n  )
        else:
            rs[n]["sig"] = -1
    rs = pd.DataFrame(rs).T
    rs = rs.fillna(0)
    rs.to_csv(fout+".txt",sep="\t")
    s = rs.loc[ts,"fc"]
    s = s.sort_values(ascending=False)
    print(s)
    with open(fout+".list","w") as fo:
        ns = [t.split("|")[0] for t in s.index]
        fo.write("\n".join(ns))
    return s.index
 

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
        print(c)
        a = cs[cs==c].index
        b = cs.index.difference(a)
        rs = checkDiff( mat[a],mat[b],"cluster_%s"%(c))
        cells.extend(a)
        gs.extend(rs)
    """
    mat = mat.loc[gs,cells]
    mat = mat.fillna(0)
    for t in mat.itertuples():
        nt = np.array(t[1:])
        nt = (nt - nt.mean())/nt.std()
        mat.loc[t[0],] = nt
    mat.to_csv("%s_clusters.txt"%pre,sep="\t")
    """


get("all")
