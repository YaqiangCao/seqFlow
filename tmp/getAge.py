import pandas as pd
import numpy as np

"""
getAge.py
2016-03-15: refer to 'The majority of primate-specific regulatory sequences are derived from transposable elements'
"""

def getDiv(f):
    mat = pd.read_table(f,index_col=None,header=0)
    nts = []
    for t in mat.itertuples():
        nt = t[4:8]
        nt = list(nt)
        nt[1] += 1
        nts.append("|".join(map(str,nt)))
    mat.index = nts
    s = mat["#milliDiv"]
    return s

def getMat(f):
    mat = pd.read_table(f,index_col=0)
    ns = mat.index
    nns = {}
    for n in ns:
        nt = n.split("|")
        nt = "|".join(nt[:4])
        nns[nt] = n
    nns = pd.Series(nns)
    return nns,mat

def preAge(s,ns,mat):
    ts = s.index.intersection(ns.index)
    s = s[ts]
    s.index = ns[ts]
    mat = mat.loc[s.index,:]
    #used for UCSC result as milliDiv, have been already multipy 1000
    mat["div"] = s / 1000.0
    return mat
    
def getAge(mat,p=2.2e-9):
    s = mat["div"]
    #simplest method as divid substituation rate
    #s = s / p / 1e6
    #Jukes-Cantor method
    s = -3.0/4*np.log(1-4.0/3*s) / p /1e6
    #mat["age"] = s
    s.to_csv("age.txt",sep="\t")


s = getDiv("UCSC_RMSK.txt")
ns,mat = getMat("hg38Rep.txt")
mat = preAge(s,ns,mat)
getAge(mat)
