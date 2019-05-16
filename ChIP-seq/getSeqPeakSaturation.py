#--coding:utf-8--
"""
getSeqPeakSaturation.py
Due to MACS require python2.7, so the programme should be run with python2.7.
"""

#sys
import os, gzip, random
from glob import glob
from datetime import datetime

#3rd
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

#global settings
import matplotlib as mpl
mpl.use("pdf")
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["figure.figsize"] = (4, 2.75)
mpl.rcParams["font.size"] = 8.0
#import seaborn as sns
import brewer2mpl
import pylab
#sns.set_style("white")
colors = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors


def callSys(cmds):
    """
    Call systematic commands without return.
    """
    for c in cmds:
        print(c)
        try:
            os.system(c)
        except:
            print("ERROR!", c)


def callPeaks(bed, toRoot="./samplingPeaks/"):
    sample = bed.split("/")[-1]
    fo = toRoot + sample
    doMACS = "macs2 callpeak -t {bed} -n {sample} -g mm --nomodel --shift 75 --extsize 150 --keep-dup 1 ".format(
        bed=bed, sample=fo)
    callSys([doMACS])
    return fo + "_peaks.narrowPeak"


def callSamplingPeaks(f, ratio, repeats=5, toRoot="./samplingBeds/",
                      pre="tmp"):
    ds = gzip.open(f, 'rt').read().split("\n")
    n = len(ds)
    nn = int(n * ratio)
    for i in range(repeats):
        random.shuffle(ds)
        nds = ds[:nn + 1]
        fout = toRoot + pre + "_%s" % ratio + "_rep%s" % i
        with open(fout, "w") as fo:
            fo.write("\n".join(nds))
        callPeaks(fout)
        callSys(["rm %s" % fout])


def getN(f):
    c = open(f).read().split("\n")
    c = [t for t in c if t.strip() != ""]
    c = len(c)
    return c


def summaryOverlaps(ref,pre="test"):
    jis = {}
    rs = {}
    fs = glob("./samplingPeaks/*.narrowPeak")
    fs.sort()
    ci = getN(ref)
    for f in fs:
        n = f.split("/")[-1].split("_")
        ratio = float("%.2f" % float(n[1]))
        if ratio not in jis:
            jis[ratio] = {}
        if ratio not in rs:
            rs[ratio] = {}
        rep = n[2]
        cj = getN(f)
        cmd = "intersectBed -a %s -b %s -u > ./tmp.bed" % (ref, f)
        callSys([cmd])
        cij = getN("./tmp.bed")
        ji = cij / 1.0 / (ci + cj - cij)
        r = cij / 1.0 / ci
        jis[ratio][rep] = ji
        rs[ratio][rep] = r
    jis = pd.DataFrame(jis)
    jis.to_csv("%s_JI.txt"%pre, sep="\t", index_label="sample")
    rs = pd.DataFrame(rs)
    rs.to_csv("%s_ratios.txt"%pre, sep="\t", index_label="sample")


def plotSummary(f):
    mat = pd.read_csv(f, index_col=0, sep="\t")
    fig, ax = pylab.subplots()
    x = list(mat.columns)
    y = list(mat.mean())
    yerr = list(mat.std())
    yerr = np.array([[t] for t in yerr])
    ax.errorbar(x, y, yerr=yerr, fmt='--o')
    ax.set_xlabel("reads re-sampling ratio")
    ax.set_ylabel(f.split("/")[-1].split("_")[-1].split(".")[0])
    pylab.savefig(f.replace(".txt", ".pdf"))


def main():
    for peakf in glob("../2.1.peaks/*.narrowPeak"):
        pre = peakf.split("/")[-1].split(".narrowPeak")[0]
        print(pre)
        bed = "../1.beds/" + pre + ".bed.gz"
        if not os.path.isfile(bed):
            continue
        for ratio in np.arange(0.1,1.0,0.1):
            callSamplingPeaks(bed,ratio)
        summaryOverlaps(peakf,"%s_stat"%pre)
        callSys(["rm samplingBeds/*","rm samplingPeaks/*"])
    for f in glob("*.txt"):
        plotSummary(f)


if __name__ == "__main__":
    startTime = datetime.now()
    main()
    elapsed = datetime.now() - startTime
    print("The process is done.")
    print("Time used:", elapsed)
