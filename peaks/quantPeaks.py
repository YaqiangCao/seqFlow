"""
quantPeaks.py
"""

import os, gzip
from glob import glob

#3rd
import pylab
import HTSeq
import brewer2mpl
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm
from joblib import Parallel, delayed

#global settings
import matplotlib as mpl
mpl.use("pdf")
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["figure.figsize"] = (4, 2.75)
mpl.rcParams["font.size"] = 10.0
sns.set_style("white")
colors = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors


def buildCovModel(readF):
    """
    readF: bed.gz
    """
    model = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    for i, line in tqdm(enumerate(gzip.open(readF, 'rt'))):
        line = line.split("\n")[0].split("\t")
        iv = HTSeq.GenomicInterval(line[0], int(line[1]), int(line[2]))
        model[iv] = i
    return model, i


def quantify(readF, peakF, fnOut):
    if os.path.isfile(fnOut):
        return  #has been generated
    print("builidng coverage model for counting")
    covModel, t = buildCovModel(readF)
    r = set()
    ds = {}
    print("counting reads in peaks")
    for line in tqdm(list(open(peakF))):
        line = line.split("\n")[0].split("\t")
        iv = HTSeq.GenomicInterval(line[0], int(line[1]), int(line[2]))
        c = set()
        for ivb, vb in covModel[iv].steps():
            try:
                c.add(vb)
            except:
                continue
        r.update(list(c))
        c = list(c)
        ds["|".join(line[:3])] = {
            "count": len(c),
            "RPKM": len(c) / 1.0 / iv.length / t * 10**9,
            "TPM": len(c)/1.0 /iv.length * 10**3,
            "length": iv.length
        }
    ds = pd.DataFrame(ds).T
    ds["TPM"] = ds["TPM"] / ds["TPM"].sum() * 10**6
    ds.to_csv(fnOut,sep="\t",index_label="peakId")


def getAllCounts():
    peakf = "../1.call/Trac_SOD_borders.bed"
    beds = glob("/mnt/data/caoy7/Projects/5.4DN_Trac/0.Pub/0.ENCODE/5.ChIP-seq/1.GM12878/3.beds/*.bed.gz")
    ds = Parallel(n_jobs=20)(delayed(quantify)( bed,peakf, "./peaksQuant/%s_peaksQuant.txt"%bed.split("/")[-1].replace(".bed.gz","") ) for bed in beds)
  


def summaryCounts():
    fs = glob("*%s*.txt"%p)
    fs.sort()
    fs.reverse()
    countDs = {}
    RPKMDs = {}
    TPMDs = {}
    for f in fs:
        mat = pd.read_csv(f,index_col=0,sep="\t")
        n = f.split("/")[-1].replace("_peaksQuant.txt","")
        countDs[n] = mat["count"]
        RPKMDs[n] = mat["RPKM"]
        TPMDs[n] = mat["TPM"]
    countDs = pd.DataFrame( countDs )
    RPKMDs = pd.DataFrame( RPKMDs )
    TPMDs = pd.DataFrame( TPMDs )
    countDs = countDs.fillna(0.0)
    RPKMDs = RPKMDs.fillna(0.0)
    TPMDs = TPMDs.fillna(0.0)
    countDs.to_csv( p + "_count.txt",sep="\t",index_label="peakId")
    RPKMDs.to_csv( p + "_RPKM.txt",sep="\t",index_label="peakId")
    TPMDs.to_csv( p + "_TPM.txt",sep="\t",index_label="peakId")


def main():
    #getAllCountsRPKM()
    summaryCountsRPKM("borders")


if __name__ == "__main__":
    main()
