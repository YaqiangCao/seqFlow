#!/usr/bin/env python2.7
#--coding:utf-8--


"""
"""



__author__="CAO Yaqiang"
__date__="2014-06-19"
__modified__=""
__email__="caoyaqiang0410@gmail.com"


#systematic library
import glob,os
import pandas as pd
import bigfloat
import numpy as np



def callSys(cmds):
    """
    Call systematic commands without return.
    """
    for c in cmds:
         print c
         try:
             os.system(c)
         except:
             print "ERROR for %s" % c


def runCentrimo():
    motifs = "/picb/molsysbio/usr/caoyaqiang/4.ENV/2.Bios/meme_4.10.0/new_motif_databases/motif_databases/HUMAN/HOCOMOCOv10_HUMAN_mono_meme_format.meme"
    fs = glob.glob("../../4.1.Fas/*.fa")
    print fs
    for fa in fs:
        fn = fa.split("/")[-1].split(".")[0]
#        if os.path.isdir(fn):
#            continue
        cmd = "centrimo --verbosity 1 -oc {out} {fa} {motifs} &".format(out=fn,fa=fa,motifs=motifs)
        callSys([cmd])


def norm(ts):
    vs = ts.values
    vs = [ 0.0-float( bigfloat.log10(bigfloat.BigFloat(v)) ) for v in vs]
    ts = pd.Series(vs,index=ts.index)
    return ts


def getMat(cut=1e-2):
    cut = 0.0 - np.log10(cut)
    #metaf = "/picb/molsysbio/usr/caoyaqiang/4.ENV/2.Bios/meme_4.10.0/new_motif_databases/motif_databases/HUMAN/HUMAN_mono_motifs.tsv"
    #metaf = "/picb/molsysbio/usr/caoyaqiang/4.ENV/2.Bios/meme_4.10.0/new_motif_databases/motif_databases/MOUSE/MOUSE_mono_motifs.tsv"
    metaf = "/picb/molsysbio/usr/caoyaqiang/4.ENV/2.Bios/meme_4.10.0/new_motif_databases/motif_databases/HUMAN/motif.txt"
    #ms = pd.read_table(metaf,index_col=0)["Transcription factor"]
    ms = pd.Series.from_csv(metaf,sep="\t")
    fs = glob.glob("*/centrimo.txt")
    ds = {}
    for f in fs:
        fn = f.split("/")[-2]
        mat = pd.read_table(f,index_col=1,skiprows=1)
        mat.index = [i.split()[0].strip() for i in mat.index] 
        s = mat[" E-value"]
        rs = [t for t in s.index if "_" in t if not t.endswith(".S")]
        s = s[rs]
        s.index = [t.split("_")[0] for t in s.index]
        rs = s.index.intersection(ms.index)
        s = s[rs]
        s.index = ms[s.index]
        ds[fn] = s
    ds = pd.DataFrame(ds)
    ds = ds.fillna(1.0)
    for c in ds.columns:
        ds[c] = norm(ds[c])
    t = ds.max(axis=1)
    t = t[t>cut]
    ds = ds.loc[t.index,:]
    ds[ds<0] = 0
    rs = [t[0] for t in ds.itertuples() if np.sum(t[1:])==0 ]
    ds = ds.drop(rs)
    ds.to_csv("e_value.txt",sep="\t",index_label="TF")


def getMat2(cut=10):
    metaf = "/picb/molsysbio/usr/caoyaqiang/4.ENV/2.Bios/meme_4.10.0/new_motif_databases/motif_databases/HUMAN/motif.txt"
    ms = pd.Series.from_csv(metaf,sep="\t")
    fs = glob.glob("*/centrimo.txt")
    ds = {}
    for f in fs:
        fn = f.split("/")[-2]
        mat = pd.read_table(f,index_col=1,skiprows=1)
        mat.index = [i.split()[0].strip() for i in mat.index] 
        s =0.0- mat["log_adj_p-value"]
        rs = [t for t in s.index if "_" in t if not t.endswith(".S")]
        s = s[rs]
        s.index = [t.split("_")[0] for t in s.index]
        rs = s.index.intersection(ms.index)
        s = s[rs]
        s.index = ms[s.index]
        ds[fn] = s
    ds = pd.DataFrame(ds)
    ds = ds.fillna(0.0)
    t = ds.min(axis=1)
    t = t[t>cut]
    ds = ds.loc[t.index,:]
    rs = [t[0] for t in ds.itertuples() if np.sum(t[1:])==0 ]
    ds = ds.drop(rs)
    ds.to_csv("e_value.txt",sep="\t",index_label="TF")



def main(  ):
    runCentrimo()
    #getMat()
    #getMat2()



main(  )
