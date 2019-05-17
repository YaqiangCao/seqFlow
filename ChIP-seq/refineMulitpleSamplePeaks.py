#!/Users/caoyaqiang/anaconda3/bin/python3.6

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

#global settings
import matplotlib as mpl
mpl.use("pdf")
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["figure.figsize"] = (4, 2.75)
mpl.rcParams["font.size"] = 10.0
sns.set_style("white")
colors = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors


def sortPeaks():
    fs = glob("../../2.1.peaks/*narrowPeak")
    for f in fs:
        n = f.split("/")[-1]
        cmd = "sortBed -i %s > ./tmp2/%s" % (f, n)
        os.system(cmd)


def getMergedSets(pattern="_-6_CD4_",
                  skips=[
                      "../3.peaks/WT_6_Th2_Neg_Rep1.narrowPeak",
                      "../3.peaks/WT_6_Th2_Neg_Rep2.narrowPeak"
                  ]):
    fs = glob("./tmp2/*%s*" % pattern)
    fs.sort()
    fs = [f for f in fs if f not in skips]
    names = [f.split("/")[-1].split(".")[0] for f in fs]
    #cmd = "multiIntersectBed -cluster -header -i %s -names %s > %s.txt "%(" ".join(fs)," ".join(names),pattern)
    cmd = "/Users/caoyaqiang/Projects/2.AID_GangRen/0.packages/bedtools2/bin/multiIntersectBed -header -i %s -names %s > %s.txt " % (
        " ".join(fs), " ".join(names), pattern)  # too big with meanless result
    print(cmd)
    os.system(cmd)


def parseMergedClusters(f, cut=1):
    mat = pd.read_csv(f, index_col=None, header=0, sep="\t")
    #select high quality peaks
    s = mat["num"]
    s = s[s > cut]
    mat = mat.loc[s.index, ]
    nmat = {}
    columns = mat.columns[5:]
    peaks = []
    for t in mat.itertuples():
        ni = "|".join(list(map(str, t[1:4])))
        s = pd.Series(t[6:], index=columns)
        nmat[ni] = s
    nmat = pd.DataFrame(nmat).T
    nmat.to_csv(f.replace(".txt", "_filter.txt"), sep="\t", index_label="peak")


def txt2bed(f):
    with open(f.replace(".txt", ".bed"), "w") as fo:
        for line in open(f).read().split("\n")[1:]:
            nline = line.split("\t")[0]
            line = nline.split("|")
            line.append(nline)
            fo.write("\t".join(line) + "\n")


def mergeBed(f):
    fo = f.replace(".bed", "_merged.bed")
    cmd = "mergeBed -d 150 -i %s > %s" % (f, fo)
    print(cmd)
    os.system(cmd)
    #select peaks, those shorter than 50 were removed
    ds = []
    for line in open(fo):
        nline = line.split("\n")[0].split("\t")
        if int(nline[2]) - int(nline[1]) > 50:
            ds.append(line)
    with open(fo, "w") as fo:
        fo.write("".join(ds))


def main():
    #sortPeaks()
    #for p in ["_-6_CD4_","_6_Th1_","_6_Th2_","_24_Th1_","_24_Th2_"]:
    #    getMergedSets( p )
    #    parseMergedClusters("./%s.txt"%p,cut=1)
    #for f in glob("*filter.txt"):
    #    txt2bed(f)
    #for f in glob("*.bed"):
    #    if "merged" not in f:
    #        mergeBed(f)
    #all sites
    getMergedSets("_")
    parseMergedClusters("./_.txt", cut=1)
    txt2bed("_.txt")
    mergeBed("_.bed")


if __name__ == "__main__":
    main()
