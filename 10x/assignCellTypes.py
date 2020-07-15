#sys
import random
from glob import glob

#3rd
import numpy as np
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed

global mat
mat = None


def assignCellType(cid, cellMarkers, perm=1000):
    """
    Assign single cell to a cell type by reference markers through permutation.
    @param cid: str, column/cell id for a matrix 
    @param cellMarkers: dict, key is cell type and values are gene id
    @perm: int, permutation times 
    @return: dict 
    """
    global mat
    s = mat[cid]
    #foreground
    ds = {}
    for c, gs in cellMarkers.items():
        ns = s[gs]
        ens = ns[ns > 0]
        ds[c] = {
            "meanExp": ns.mean(),
            "expressedMeanExp": ens.mean(),
            "expressedGenes": len(ens) / len(ns),
            "markerGenes": len(gs)
        }
    bgs = {}
    #permutation background
    #for i in tqdm(range(perm)):
    for i in range(perm):
        v = s.values
        np.random.shuffle(v)
        s = pd.Series(v, index=s.index)
        for c, gs in cellMarkers.items():
            ns = s[gs]
            ens = ns[ns > 0]
            if c not in bgs:
                bgs[c] = {
                    "meanExp": [],
                    "expressedMeanExp": [],
                    "expressedGenes": []
                }
            bgs[c]["meanExp"].append(ns.mean())
            bgs[c]["expressedMeanExp"].append(ens.mean())
            bgs[c]["expressedGenes"].append(len(ens) / len(ns))
    #caculation the fdr and enrichment score for each cell type
    for c in bgs.keys():
        sa = np.array(bgs[c]["meanExp"])
        fdr = len(sa[sa > ds[c]["meanExp"]]) / len(sa)
        ds[c]["fdr"] = fdr
        ds[c]["ES"] = ds[c]["meanExp"] / np.mean(bgs[c]["meanExp"])
        ds[c]["bg_meanExp"] = np.mean(bgs[c]["meanExp"])
        ds[c]["bg_expressedMeanExp"] = np.mean(bgs[c]["expressedMeanExp"])
        ds[c]["bg_expressedGenes"] = np.mean(bgs[c]["expressedGenes"])
    ds = pd.DataFrame(ds).T
    ds = ds.fillna(0.0)
    ds.to_csv("1.initCellTypes/%s.txt" % cid, sep="\t")


def get(gmtf, expf):
    global mat
    cs = {}
    for line in open(gmtf):
        line = line.split("\n")[0].split("\t")
        cs[line[0]] = line[2:]
    mat = pd.read_csv(expf, sep="\t", index_col=0)
    mat.index = [t.split("|")[0] for t in mat.index]
    #others genes
    ns = []
    for k, v in cs.items():
        ns.extend(v)
    ns = mat.index.difference(ns)
    ns = list(ns)
    random.shuffle(ns)
    for k,v in cs.items():
        v = list(mat.index.intersection(v))
        cs[k] = v
    cs["random_others"] = ns[:500]
    Parallel(n_jobs=60,
             backend="multiprocessing")(delayed(assignCellType)(c, cs)
                                        for c in tqdm(mat.columns))


def summary(dir="1.initCellTypes/",fdr=0.005,es=2,er=0.2):
    fs = glob(dir+"/*.txt")
    ds = {}
    for f in tqdm(fs):
        cid = f.split("/")[-1].split(".txt")[0]
        mat = pd.read_csv(f,index_col=0,sep="\t")
        a = mat["fdr"]
        a = a[a<=fdr]
        if len(a) == 0:
            ds[cid] = {"potentialCellNumbers":0}
            continue
        b = mat.loc[a.index,"ES"]
        b = b[b>=es]
        if len(b) == 0:
            ds[cid] = {"potentialCellNumbers":0}
            continue
        c = mat.loc[b.index,"expressedGenes"]
        c = c[c>=er]
        if len(c) == 0:
            ds[cid] = {"potentialCellNumbers":0}
            continue
        d = mat.loc[c.index,"meanExp"]
        d = d/mat.loc["random_others","meanExp"]
        d = d[d>=es]
        if len(d) == 0:
            ds[cid] = {"potentialCellNumbers":0}
            continue
        e = mat.loc[d.index,"ES"]
        e = e.sort_values(inplace=False,ascending=False)
        cts = ",".join(list(e.index))
        ds[cid] = {
            "potentialCellNumbers": len(e),
            "potentialCellTypes":",".join(list(e.index)),
            "ES":",".join(list(map(str,mat.loc[e.index,"ES"]))),
            "meanExp":",".join(list(map(str,mat.loc[e.index,"meanExp"]))),
            "expressedGenes":",".join(list(map(str,mat.loc[e.index,"expressedGenes"]))),
            }
    ds = pd.DataFrame(ds).T
    ds = ds.fillna(0)
    ds.to_csv("1_CellTypes.txt",sep="\t") 


#gmtf = "../0.ref/GSE112393_MouseTestis_cellMarkers.gmt"
#expf = "../6.mergedAll/all_normed.txt"
#get(gmtf, expf)
summary()
