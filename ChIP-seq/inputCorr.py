#!/usr/bin/env python2.7
#--coding:utf-8--
"""
inputCorr.py
Estimate the normalization factor and then normalize the ChIP data.
"""

__author__ = "CAO Yaqiang"
__date__ = ""
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import os, glob, gzip, sys
from datetime import datetime
from collections import Counter

#plotting settings
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

#3rd library
import HTSeq
import numpy as np
import pandas as pd
from joblib import Parallel, delayed


def commandFlush(r):
    """
    One line flush.
    """
    sys.stdout.write("\r%s" % r)
    sys.stdout.flush()


def getRegions(chrom1, chrom2):
    """
    Get overlapped chromosomes and the start and end regions.
    """
    chroms = {}
    cs = set(chrom1.keys()).intersection(chrom2.keys())
    cs = list(cs)
    for c in cs:
        s = min(chrom1[c]["s"], chrom2[c]["s"])
        e = max(chrom1[c]["e"], chrom2[c]["e"])
        chroms[c] = {"s": s, "e": e}
    return chroms


def getCount(model, iv, r):
    """
    Get reads count from model in defined region.
    """
    c = 0
    for ivb, value in model[iv].steps():
        c += ivb.length * value
    #get upper bound interger
    c = np.ceil(float(c) / r)
    return c


def bed2GModel(bed, ext=150):
    """
    Bed format, gzip or not into HTSeq.GenomicArray
    """
    if bed.endswith(".gz"):
        f = gzip.open(bed, "rb")
    else:
        f = open(bed)
    print datetime.now(), "Start building model for %s" % bed
    model = HTSeq.GenomicArray("auto", stranded=False)
    chroms = {}
    for i, line in enumerate(f):
        if i % 10000 == 0:
            report = "%s lines genome signal read." % i
            commandFlush(report)
        line = line.split("\n")[0].split("\t")
        try:
            chrom = line[0]
            s = int(line[1])
            e = int(line[2])
            strand = line[5]
        except:
            continue
        if s - e < ext:
            if strand == "+":
                e = s + ext
            else:
                s = e - ext
                if s < 0:
                    s = 0
        iv = HTSeq.GenomicInterval(chrom, s, e, strand)
        model[iv] += 1
        if chrom not in chroms:
            chroms[chrom] = {"s": s, "e": e}
        if s < chroms[chrom]["s"]:
            chroms[chrom]["s"] = s
        if e > chroms[chrom]["e"]:
            chroms[chrom]["e"] = e

    print
    print datetime.now(), "Model built for %s" % bed
    return model, chroms, i, bed.split("/")[-1].split(".bed.gz")[0]


def binData(tmodel, cmodel, tr, cr, chroms, binsize):
    """
    Get bins for the genome and get the reads count.
    tr is read length for treatment
    cr is read length for control
    """
    ts = []
    cs = []
    for chrom in chroms.keys():
        s = chroms[chrom]["s"]
        e = chroms[chrom]["e"]
        bins = (e - s) / binsize
        for i in xrange(1, bins):
            iv = HTSeq.GenomicInterval(chrom, s + binsize * (i - 1),
                                       s + binsize * i)
            countT = getCount(tmodel, iv, tr)
            countC = getCount(cmodel, iv, cr)
            if countT + countC == 0:
                #if countT * countC == 0:
                continue
            else:
                ts.append(countT)
                cs.append(countC)
    ts, cs = np.array(ts), np.array(cs)
    return ts, cs


def estNorm(tmodel, cmodel, tname, cname, tr, cr, chroms, binvecs, quant=0.85):
    """
    Estimate the normalization factor
    """
    print datetime.now(), "Start estimating the normalization facotr."
    brs = []
    counts = []
    est = None
    for binsize in binvecs:
        #marks this window size found the best reads count cutoff
        subflag = 0
        ts, cs = binData(tmodel, cmodel, tr, cr, chroms, binsize)
        total = ts + cs
        ns = pd.Series(Counter(total))
        ns = np.cumsum(ns) / total.shape[0]
        clow = int(max(1, max(ns[ns <= quant].index)))
        chigh = int(max(np.max(ts), np.max(cs))) + 1
        trs = []
        for t in xrange(clow, chigh):
            #report = "Estimating the normalization facotr for bin size of %s and count cutoff of %s" % (
            #    binsize, t
            #)
            #commandFlush(report)
            ind = np.where(total < t)[0]
            if len(ind) == 0:
                continue
            tmpts, tmpcs = np.sum(ts[ind]), np.sum(cs[ind])
            if tmpts * tmpcs == 0:
                continue
            r = float(tmpts) / tmpcs
            report = "Estimating the normalization facotr for {tname} vs {cname}: binSize={binSize},countCutoff={countCutoff},chipRatio={chipRatio}, inputRatio={inputRatio}, estR={estR}".format(
                tname=tname,
                cname=cname,
                binSize=binsize,
                countCutoff=t,
                chipRatio=tmpts / np.sum(ts) * 100,
                inputRatio=tmpcs / np.sum(cs) * 100,
                estR=r)
            print report
            if len(trs) > 0 and r > trs[-1]:
                subflag = 1
                break
            trs.append(r)
        if subflag == 1:
            est = r
            if len(brs) > 0 and r > brs[-1]:
                break
            brs.append(r)
        print
    print
    #modified ts and cs for nice plot
    #total = ts + cs
    #ns = pd.Series(Counter(total))
    #ns = np.cumsum(ns) / total.shape[0]
    #chigh = max(ns[ns <= 0.95].index)
    #ind = np.where(total <= chigh)[0]
    ts = ts[ind]
    cs = cs[ind]
    return binsize, t, ts, cs, est


def getNCIS(t, c, fr=150, pre=None, board=True, quant=0.85):
    """
    t = [tmodel, tchroms, tn, tr,tname] 
    c = [cmodel, cchroms, cn, cr,cname] 
    Multiple processing.
    """
    binVecs = []
    if board:
        binVecs = [5000, 10000, 20000, 50000]
    else:
        binVecs = [100, 200, 500, 1000, 2000, 5000, 10000]
    tmodel, tchroms, tn, tname = t[0], t[1], t[2], t[3]
    cmodel, cchroms, cn, cname = c[0], c[1], c[2], c[3]
    tr = fr
    cr = fr
    depthRatio = float(tn) / cn
    chroms = getRegions(tchroms, cchroms)
    binSize, usedCount, ts, cs, est = estNorm(tmodel, cmodel, tname, cname, tr,
                                              cr, chroms, binVecs, quant)
    if est == None:
        print "NO BINDING SIGNAL FOR %s and %s!!" % (tname, cname)
        return
    PI = est / depthRatio
    #print binSize, usedCount, est, PI, depthRatio
    #plot scatter plot and estimation of PI and
    fig, ax = pylab.subplots()
    df = pd.DataFrame({"Background": cs, "ChIP": ts})
    df.plot(
        kind="hexbin",
        x="Background",
        y="ChIP",
        ax=ax,
        colormap=pylab.cm.BuPu,
        title="PI:%.3f,binSize:%s\nusedCountCufoff:%d" %
        (PI, binSize, usedCount),  #reduce_C_function=np.median,
        reduce_C_function=np.max,
        bins="log")
    ax.plot([0, np.max(cs)], [0, depthRatio * np.max(cs)],
            label="DR:%.3f" % depthRatio,
            c="k")
    ax.plot([0, np.max(cs)], [0, est * np.max(cs)],
            label="estR:%.3f" % est,
            c="b")
    ax.legend(loc="best")
    pylab.savefig(tname + "_est.pdf")
    return tname, binSize, usedCount, est, PI, depthRatio


def getInputCorr(tmodel, cmodel, tname, cname, estR):
    """
    Normalize the signal to input as s = max(0,s-R*i )
    """
    print datetime.now(), "start %s normalize to %s" % (tname, cname)
    for i, (iv, s) in enumerate(tmodel.steps()):
        if i % 10000 == 0:
            report = "%s records of %s normalized to %s" % (i, tname, cname)
            commandFlush(report)
        if s == 0:
            continue
        s = max(0, s - estR * np.array(list(cmodel[iv])).mean())
        tmodel[iv] = s
    print
    print datetime.now(), "%s finished normalize to %s" % (tname, cname)
    return tmodel


def writeBg(bgz, model):
    f = gzip.open(bgz, "wb")
    for i, (iv, s) in enumerate(model.steps()):
        if i % 10000 == 0:
            report = "%s records writen for %s" % (i, bgz)
            commandFlush(report)
        if s == 0:
            continue
        line = [iv.chrom, iv.start, iv.end, s]
        f.write("\t".join(map(str, line)) + "\n")
    f.close()
    print
    print datetime.now(), "writing finished for %s" % bgz


def preDs(beds):
    """
    Prepare the datasets as chip-input. 
    """
    ds = {}
    for bed in beds:
        fn = bed.split("/")[-1].split(".bed.gz")[0].split("_")
        c, t = fn[0], fn[-1]
        if c not in ds:
            ds[c] = {"chip": [], "input": None}
        if t == "DNase":
            continue
        if t == "Input":
            ds[c]["input"] = bed
        else:
            ds[c]["chip"].append(bed)
    for c in ds.keys():
        if ds[c]["input"] == None:
            del ds[c]
    return ds


def minputCorr(input, chips, pre, fr=150, board=True, quant=0.9):
    fout = "%s_estR.txt" % pre
    if os.path.exists(fout):
        print datetime.now(), "%s has been generated, return" % fout
        return
    print datetime.now(), "%s started!" % pre
    data = {}
    c = bed2GModel(input, fr)
    for chip in chips:
        fout2 = chip.split("/")[-1].split(".bed.gz")[0] + "_n.bg.gz"
        if os.path.exists(fout2):
            print "%s has been generated, return" % fout2
            continue
        t = bed2GModel(chip, fr)
        tname, binSize, usedCount, estR, PI, depthRatio = getNCIS(t,
                                                                  c,
                                                                  board=board,
                                                                  quant=quant)
        tmodel = getInputCorr(t[0], c[0], t[-1], c[-1], estR)
        writeBg(fout2, tmodel)
        data[tname] = {
            "binSize": binSize,
            "usedCountCutoff": usedCount,
            "estR": estR,
            "PI": PI,
            "DepthRatio": depthRatio
        }
    data = pd.DataFrame(data).T
    data.to_csv(fout, sep="\t")
    print datetime.now(), "%s finished!" % pre


def main():
    beds = glob.glob(
        "../../../1.ProcessedTagAlign/2.BEDhg38/1.Consolidated/*.bed.gz")
    ds = preDs(beds)
    Parallel(n_jobs=25)(
        delayed(minputCorr)(ds[cell]["input"], ds[cell]["chip"], cell)
        for cell in ds.keys())


if __name__ == '__main__':
    start_time = datetime.now()
    main()
    elapsed = datetime.now() - start_time
    fn = os.path.basename(__file__)
    print datetime.now(), "The process is done for %s,time used:%s" % (fn,
                                                                       elapsed)
