#!/usr/bin/env python2.7
#--coding:utf-8--
"""
chiapetPerRep.py
2016-12-12: permutation for genomic interactions and then caculating the confidence and support.
"""

#sys
import os, random
from glob import glob
from datetime import datetime

#3rd library
import HTSeq
import numpy as np
import matplotlib as mpl
mpl.use( "pdf" )
import seaborn as sns
mpl.rcParams[ "pdf.fonttype" ] = 42
mpl.rcParams[ "figure.figsize" ] = ( 4,2.75 )
mpl.rcParams[ "figure.dpi" ] = 100
mpl.rcParams[ "savefig.transparent" ] = True
mpl.rcParams[ "savefig.bbox" ] = "tight"
mpl.rcParams[ "font.size" ] = 10.0
mpl.rcParams[ "font.sans-serif" ] = "Arial"
mpl.rcParams[ "savefig.format" ] = "pdf"
import pylab
sns.set_style( "white" )
import numpy as np
import pandas as pd
from joblib import Parallel,delayed
import brewer2mpl
colors =  brewer2mpl.get_map('Set1', 'qualitative', 9).mpl_colors    
colors.extend(brewer2mpl.get_map('Set3', 'qualitative', 8).mpl_colors)


#global values
global repModel
global gs


def readRep(f):
    repModel = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    for i, line in enumerate(open(f)):
        if i == 0:
            continue
        line = line.split("\n")[0].split("\t")
        locus = [line[0], line[1], line[2], line[3]]
        try:
            iv = HTSeq.GenomicInterval(locus[1], int(locus[2]), int(locus[3]))
            repModel[iv] += locus[0] + "|" + line[6]
        except:
            continue
    return repModel


def getGenomeSize():
    f = "/picb/molsysbio/usr/caoyaqiang/1.Projects/17.Alu/1.GenomeReference/1.hg38/1.hg38_Sequence/hg38.chrom.sizes.mango"
    gs = pd.Series.from_csv(f, sep="\t")
    return gs


def getBedpe(f):
    fb = f.split("/")[-1].split(".")[0] + ".bedpe"
    with open(fb, "w") as fo:
        for line in open(f):
            line = line.split("\n")[0].split("\t")
            nline = line[:7]
            fo.write("\t".join(nline) + "\n")
    return fb


def getRandomRegion(chrom, start, end):
    global gs
    start, end = int(start), int(end)
    nchrom = gs.index[random.randint(0, len(gs) - 1)]
    nstart = random.randint(0, gs[nchrom] - (end - start))
    nend = nstart + (end - start)
    return [nchrom, nstart, nend]


def randomBedpe(bedpe):
    tmp = str(random.random())
    with open(tmp, "w") as f:
        for line in open(bedpe):
            try:
                line = line.split("\n")[0].split("\t")
                nline = []
                a = getRandomRegion(line[0], line[1], line[2])
                b = getRandomRegion(line[3], line[4], line[5])
                nline.extend(a)
                nline.extend(b)
                f.write("\t".join(map(str, nline)) + "\n")
            except:
                continue
    return tmp


def getRepOverlap(iv):
    global repModel
    reps = set()
    for niv, value in list(repModel[iv].steps()):
        reps.update(value)
    return reps


def mapBedpeToRep(bedpe):
    ds = []
    for i, line in enumerate(open(bedpe)):
        if i == 0:
            continue
        line = line.split("\n")[0].split("\t")
        a = HTSeq.GenomicInterval(line[0], int(line[1]), int(line[2]))
        b = HTSeq.GenomicInterval(line[3], int(line[4]), int(line[5]))
        ar = getRepOverlap(a)
        br = getRepOverlap(b)
        if len(ar) == 0 and len(br) == 0:
            continue
        ds.append([ar, br])
    return ds


def preRepPairs(ds):
    rs = []
    for d in ds:
        a, b = d[0], d[1]
        a = [t.split("|")[-1] for t in a]
        b = [t.split("|")[-1] for t in b]
        a, b = list(set(a)), list(set(b))
        if len(a) == 0 and len(b) == 0:
            continue
        rs.append([a, b])
    return rs


def cusapir(rs, term=("L2", "MIR")):
    #only caculating for ["MIR","L2"]
    ds = {}
    ca, cb = 0, 0
    a, b = term[0], term[1]
    for r in rs:
        if (a in r[0] and b in r[1]) or (a in r[1] and b in r[0]):
            ca += 1
        if (a in r[0]) or (a in r[1]):
            cb += 1
    ca, cb = float(ca) / len(rs), float(cb) / len(rs)
    #support
    s = ca
    #confidence
    c = ca / cb
    return s, c


def sPerTest(i, fb):
    print "%s permutation for %s" % (i, fb)
    tmpf = randomBedpe(fb)
    ds = mapBedpeToRep(tmpf)
    rs = preRepPairs(ds)
    s, c = cusapir(rs)
    os.remove(tmpf)
    print "%s permutation for %s, support:%s,confidence:%s" % (i, fb, s, c)
    return s, c


def getMIRL2():
    global repModel, gs
    repF = "/picb/molsysbio/usr/caoyaqiang/1.Projects/17.Alu/1.GenomeReference/1.hg38/6.Repeats_Annotation/1.GenomicLocation/hg38Rep.txt"
    repModel = readRep(repF)
    gs = getGenomeSize()
    for f in glob("../2.map2rep/*.pgl"):
        print f
        pre = f.split("/")[-1].split(".")[0]
        fb = getBedpe(f)
        ds = Parallel(n_jobs=20)(delayed(sPerTest)(i, fb) for i in xrange(1000))
        data = {}
        for i, d in enumerate(ds):
            data[i] = {"support": d[0], "confidence": d[1]}
        data = pd.DataFrame(data).T
        data.to_csv(pre + ".txt", index_label="permutationId", sep="\t")


def plotD(data,pre,scut=0.2,ccut=0.5):
    width = 0.6
    fig, ax = pylab.subplots(figsize=(2.5,1))
    rs = list(data.index)
    rs.reverse()
    data = data.loc[rs,]
    samples = list( data.index )
    ind = np.arange( len( data.index ) )
    s = data["support"].values
    ax.bar( ind,s,width,color="gray" )
    ax.set_ylabel( "support" )
    ax.set_xticks( ind + width )
    ax.set_xticklabels( samples,rotation=90,ha="right",fontsize=5  )
    ax.axhline(y=scut,linestyle="--",color="gray")
    ax2 = ax.twinx(  )
    c = data["confidence"].values
    ax2.plot( ind+width/2,c,color=colors[1] )
    for t in ax2.get_yticklabels(  ):
        t.set_color( colors[1] )
    ax2.axhline(y=ccut,linestyle="--",color=colors[1])
    ax2.set_ylabel( "confidence" )
    pylab.savefig( "Aprori_%s.pdf"%pre )


def plotSS(target="L2->MIR",xcut=0.2,ycut=0.5):
    ds = {} 
    for f in glob("*.txt"):
        pre = f.split(".")[0]
        mat = pd.read_table(f,index_col=0)
        ds[pre] = {"support":mat["support"].mean(),"confidence":mat["confidence"].mean()} 
    ds = pd.DataFrame(ds).T
    plotD(ds,target)



if __name__ == '__main__':
    start_time = datetime.now()
    getMIRL2()
    plotSS()
    elapsed = datetime.now() - start_time
    print "The process is done"
    print "Time used:", elapsed
