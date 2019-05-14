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


def getFRiPCounts(readF, peakF, fnOut):
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
            "RPM": len(c) / 1.0 / t * 10**6,
            "macsQ": line[8]
        }
    with open(fnOut, "w") as fo:
        line = "#TotalReads\t%s\n#ReadsInPeaks\t%s\n#FRiP\t%s\n" % (
            t, len(r), len(r) / 1.0 / t)
        fo.write(line)
        line = "\t".join(["peak", "count", "RPM", "macsQ"]) + "\n"
        fo.write(line)
        for k, v in ds.items():
            line = k + "\t" + "\t".join(list(map(str, v.values()))) + "\n"
            fo.write(line)


def getAllFRiPCounts():
    for f in glob("../3.peaks/*.narrowPeak"):
        pre = f.split("/")[-1].split(".")[0]
        bed = "../1.beds/" + pre + ".bed.gz"
        if not os.path.isfile(bed):
            continue
        fnOut = "1.FRiP_counts_RPM/%s_peaksQuant.txt" % pre
        if os.path.isfile(fnOut):
            continue  #has been generated
        print(bed, f, fnOut)
        getFRiPCounts(bed, f, fnOut)


def summaryFRiP():
    fs = glob("./1.FRiP_counts_RPM/*.txt")
    ds = {}
    for f in fs:
        n = "_".join(f.split("/")[-1].split("_")[:-1])
        rs = open(f).read().split("\n")
        rs = [r for r in rs if r.startswith("#")]
        frip = rs[-1].split("\n")[0].split("\t")[1]
        ds[n] = {
            "FRiP": frip,
            "mouseType": n.split("_")[0],
            "timeStamp": n.split("_")[1],
            "cellType": n.split("_")[2],
            "GFP": n.split("_")[3],
            "rep": n.split("_")[-1]
        }
    ds = pd.DataFrame(ds).T
    s = list(ds.index)
    s.sort()
    ds = ds.loc[s]
    ds.to_csv("1_FRiP_summary.txt", sep="\t", index_label="Sample")


def plotFRiP():
    ds = pd.read_csv("1_FRiP_summary.txt", sep="\t", index_col=0)
    fig, axs = pylab.subplots(1, 3, figsize=(12, 3), sharey=True)
    #for i, time in enumerate(list(set(ds["timeStamp"]))):
    for i, time in enumerate([-6,6,24]):
        s = ds["timeStamp"]
        s = s[s == time].index
        nd = ds.loc[s, ]
        #seperate AID and WT
        sa = [c for c in nd.index if "AID" in c]
        axs[i].bar(range(len(sa)),
                   nd.loc[sa, "FRiP"],
                   width=1,
                   color=colors[0])
        sb = [c for c in nd.index if "AID" not in c]
        axs[i].bar(range(len(sa) + 1,
                         len(sa) + 1 + len(sb)),
                   nd.loc[sb, "FRiP"],
                   width=1,
                   color=colors[1])
        x = range(nd.shape[0] + 1)
        labels = sa
        labels.append("")
        labels.extend(sb)
        #general settings
        axs[i].set_title(str(time) + "h")
        axs[i].set_ylim([0, 0.18])
        #axs[i].set_xlim([0,max(x)])
        axs[i].set_xticks(x)
        axs[i].set_xticklabels(labels, rotation=45, fontsize=5, ha="right")
        sns.despine(ax=axs[i])
        if i == 0:
            axs[i].set_ylabel("FRiP")
    fig.tight_layout()
    pylab.subplots_adjust(wspace=0.05)
    pylab.savefig("1_FRiP.pdf")
    pylab.savefig("1_FRiP.png")


def main():
    #getAllFRiPCounts()
    #summaryFRiP()
    plotFRiP()


if __name__ == "__main__":
    main()
