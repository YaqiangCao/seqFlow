#!/usr/bin/env python2.7
#--coding:utf-8--
"""
tracPreBam.py
2019-05-23: basically finished
2019-05-28: updated as bedpe stats
"""

__author__ = "CAO Yaqiang"
__date__ = "2019-05-23"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import os, time, commands, gzip
from glob import glob
from datetime import datetime
from collections import Counter

#3rd library
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

#this
#from utils import getLogger, callSys, PET
from utils import *

#global settings
#logger
date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
logger = getLogger(fn=os.getcwd() + "/" + date.strip() + "_" +
                   os.path.basename(__file__) + ".log")


def thinBedpe(f):
    redus = set()
    with open(f + ".2", "w") as fo:
        for i, line in enumerate(open(f)):
            line = line.split("\n")[0].split("\t")
            #remove redudant PETs
            t = line[:6]
            t.extend([line[8], line[9]])
            t = tuple(t)
            if t in redus:
                continue
            redus.add(t)
            #remove the chr1_, chr2_dask
            if "_" in line[0] or "_" in line[3]:
                continue
            #shroten the name
            line[6] = str(i)
            fo.write("\t".join(line) + "\n")
    cmds = ["mv %s %s" % (f + ".2", f), "gzip %s" % f]
    callSys(cmds, logger)


def bam2Bedpe(bam, bedpe, mapq=10):
    fd = os.path.splitext(bedpe)[0]
    d = os.path.dirname(bedpe)
    if not os.path.exists(d):
        os.mkdir(d)
    tmpbam = fd + ".2.bam"
    rmunmaped = "samtools view -q 10 -b -F 4 {} >> {}".format(bam, tmpbam)
    callSys([rmunmaped], logger)
    bam2bedpe = "bamToBed -bedpe -i {bam} > {bedpe}".format(bam=tmpbam,
                                                            bedpe=bedpe)
    logger.info(bam2bedpe)
    status, output = commands.getstatusoutput(bam2bedpe)
    rmbam = "rm {}".format(tmpbam)
    callSys([rmbam], logger)


def convert():
    """
    Batch converting from bam to bedpe.
    """
    bams = glob("../3.mapping/*/*.bam")
    ds = []
    for bam in bams:
        bai = bam.replace(".bam", ".bai")
        if not os.path.isfile(bai):
            continue
        nb = "./" + bam.split("/")[-1].replace(".bam", ".bedpe")
        if os.path.isfile(nb):
            logger.info("%s has been generated. return." % nb)
            continue
        if os.path.isfile(nb + ".gz"):
            logger.info("%s has been generated. return." % nb + ".gz")
            continue
        ds.append([bam, nb])
    Parallel(n_jobs=10)(delayed(bam2Bedpe)(t[0], t[1]) for t in ds)
    fs = glob("*.bedpe")
    Parallel(n_jobs=10)(delayed(thinBedpe)(f) for f in fs)


def getStat(f):
    """
    Get the name, total PETs, PETs distance mean, distance std for a bedpe file.
    """
    print(f)
    ds = []
    n = f.split("/")[-1].split(".bedpe")[0]
    for i, line in enumerate(gzip.open(f)):
        line = line.split("\n")[0].split("\t")
        s = min(int(line[1]), int(line[4]))
        e = max(int(line[2]), int(line[5]))
        d = e -s
        ds.append(d)
    ds = np.array(ds)
    return n, i, ds.mean(), ds.std()


def stat():
    """
    Caculating the basic stats, like PETs number and mean distance between PET.
    """
    fs = glob("*.bedpe.gz")
    data = Parallel(n_jobs=50)(delayed(getStat)(f) for f in fs)
    ds = {}
    for d in data:
        ds[d[0]] = {"PETs": d[1], "distance mean": d[2], "distance std": d[3]}
    ds = pd.DataFrame(ds).T
    ds.to_csv("stat.txt", sep="\t")


def plotStat(cut=20000):
    """
    Plot some basic stats. Cut is the number of PETs limitation.
    """
    ds = pd.read_table("stat.txt", sep="\t", index_col=0)
    fig, ax = pylab.subplots()
    sns.boxplot(ds["distance mean"])
    pylab.savefig("stat_distance_mean.pdf")
    fig, ax = pylab.subplots()
    sns.boxplot(np.log2(ds["PETs"]))
    pylab.savefig("stat_library_complexity.pdf")
    s = ds["PETs"]
    s = s[s >= cut]
    ns = list(s.index)
    ns = ["_".join(n.split("_")[:-1]) for n in ns]
    ns = pd.Series(Counter(ns))
    rs = pd.Series(Counter(["_".join(n.split("_")[:-1]) for n in ds.index]))
    print(ns)
    print(rs)
    fig, ax = pylab.subplots(figsize=(2 * 0.8, 2.75 * 0.8))
    x = np.array(range(len(ns.index)))
    ax.bar(x, rs.values, label="raw", width=0.3, color=colors[0])
    ax.bar(x + 0.3,
           ns.values,
           label="after filtering",
           width=0.3,
           color=colors[1])
    ax.set_ylabel("Cell Number", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(list(ns.index), rotation=45, fontsize=8)
    ax.legend(fontsize=8)
    ax.set_title("Filtering requiring PETs > %s"%cut,fontsize=8)
    #fig.tight_layout()
    pylab.savefig("cell_numbers.pdf")


def main():
    #convert()
    stat()
    #plotStat()


if __name__ == '__main__':
    start_time = datetime.now()
    main()
    elapsed = datetime.now() - start_time
    print "The process is done"
    print "Time used:", elapsed
