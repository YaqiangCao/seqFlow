#!/usr/bin/env python
#--coding:utf-8--
"""
scMNaseTss.py
2019-06-12: remove redundancy added.
2019-06-27: only activate genes
"""

__author__ = "CAO Yaqiang"
__date__ = "2019-05-29"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import os, time, commands, gzip
from glob import glob
from datetime import datetime
from collections import Counter

#3rd library
import HTSeq
import numpy as np
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed

#this
from utils import *

#global settings
#logger
date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
logger = getLogger(fn=os.getcwd() + "/" + date.strip() + "_" +
                   os.path.basename(__file__) + ".log")

#TSS file
"""
example, last column is -1 means negative strand
chr1	3205901	3671498	ENSMUSG00000051951	Xkr4	-1
chr1	3999557	4409241	ENSMUSG00000025900	Rp1	-1
chr1	4490931	4497354	ENSMUSG00000025902	Sox17	-1
"""
TSS = "/home/caoy7/caoy7/Projects/0.Reference/2.mm10/2.annotations/mm10_biomart_tss_tts_proteinCoding.bed"


def buildCovModel(readF, dfilter=[80, 140, 180], mapq=1):
    """
    Building Genome Coverage profile for MNase-seq data based on HTSeq.

    Parameters
    ---
    readF: str,bedpe.gz
    dfilter: list, distance to determin conical and particle
    mapq: int, MAPQ cutoff to remove PETs.

    Returns
    ---
    non-Redundant PETs number, int
    Conical nucleosome PETs number, int
    Particle PETs number, int
    Conical nucleosome PETs coverage, HTSeq.GenomeicArray
    Particle PETs coverage, HTSeq.GenomicArray
    """
    print("building models for %s" % readF)
    n = readF.split('/')[-1].split(".bedpe.gz")[0]
    modelCn = HTSeq.GenomicArray("auto", stranded=False)
    modelSp = HTSeq.GenomicArray("auto", stranded=False)
    cn, sp = 0, 0
    reds = set()
    for i, line in enumerate(gzip.open(readF, 'rt')):
        if i % 10000 == 0:
            report = "%s lines genome signal read." % i
            cFlush(report)
        line = line.split("\n")[0].split("\t")
        if len(line) < 7:
            continue
        if line[0] != line[3]:
            continue
        if int(line[7]) < mapq:
            continue
        s = min(int(line[1]), int(line[4]))
        e = max(int(line[2]), int(line[5]))
        d = e - s
        r = (line[0], s, e)
        if r in reds:
            continue
        else:
            reds.add(r)
        m = (s + e) / 2
        iv = HTSeq.GenomicInterval(line[0], m, m + 1)
        if d <= dfilter[0]:  #sp
            sp += 1
            modelSp[iv] += 1
        if d >= dfilter[1] and d <= dfilter[2]:
            cn += 1
            modelCn[iv] += 1
    return len(reds), cn, sp, modelCn, modelSp


def getProfiles(t, model, gs=None, ext=2000, bin=1, skipZero=True):
    """
    Get profiles from the converge model for TSS.

    Parameters
    ---
    t: int, total PETs number
    model: HTSeq.GenomicArray
    gs: set, gnames to use
    ext: int, extension bp from the TSS
    bin: int, binsize
    skipZero: bool, whether to remove zero signal record

    Returns
    ---
    pd.Series with index indicate location and value indicate signal level
    """
    profile = np.zeros(2 * ext / bin)
    for line in open(TSS):
        line = line.split("\n")[0].split("\t")
        if gs is not None and line[4] not in gs:
            continue
        print(line[4])
        if line[-1] == "-1":
            m = int(line[2])
        else:
            m = int(line[1])
        s = m - ext
        if s < 0:
            continue
        e = m + ext
        iv = HTSeq.GenomicInterval(line[0], s, e, ".")
        cvg = np.fromiter(model[iv], dtype='i', count=2 * ext)
        cvg = cvg.reshape((len(cvg) / bin, bin))
        cvg = np.mean(cvg, axis=1)
        cvg = cvg / 1.0 / t * 10**6
        if skipZero and np.sum(cvg) == 0:
            continue
        if line[-1] == "-1":
            profile += cvg[::-1]
        else:
            profile += cvg
    s = pd.Series(profile, index=np.arange(0 - ext, ext, bin))
    return s


def getCnSpProfiles(f,gs=None, ext=2000, bin=1, skipZero=True, todir="data"):
    fosp = "%s/%s_sp.txt" % (todir, f.split("/")[-1].split(".bedpe")[0])
    focn = "%s/%s_cn.txt" % (todir, f.split("/")[-1].split(".bedpe")[0])
    if not os.path.exists(todir):
        os.mkdir(todir)
    if os.path.isfile(fosp) and os.path.isfile(focn):
        return
    tot, cn, sp, modelCn, modelSp = buildCovModel(f)
    cnp = getProfiles(cn, modelCn,gs, ext, bin, skipZero)
    spp = getProfiles(sp, modelSp,gs, ext, bin, skipZero)
    cnp.to_csv(focn, sep="\t")
    spp.to_csv(fosp, sep="\t")


def smooth(x, window_len=101, window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."

    if window_len < 3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    s = np.r_[x[window_len - 1:0:-1], x, x[-1:-window_len:-1]]
    if window == 'flat':  #moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.' + window + '(window_len)')

    y = np.convolve(w / w.sum(), s, mode='valid')
    return y


def plotProfileBin(fcn, fsp, fout, bin=5):
    """
    MNase-seq profile plot, left y-axis is cn, right x-axis is sp
    """
    fig, ax = pylab.subplots(figsize=(2, 2.75 * 0.8))
    ax2 = ax.twinx()
    ss = pd.Series.from_csv(fcn, sep="\t")
    x = np.arange(ss.index[0], ss.index[-1], bin)
    ss = ss.values.reshape((len(ss) / bin, bin))
    ss = np.sum(ss, axis=1)
    ax.plot(x, ss, color=colors[0], linewidth=1)
    ax.set_ylabel("Nucleosome density")
    ax.set_xlabel("Distance from TSS")
    ax.set_title(fout)
    for t in ax.get_yticklabels():
        t.set_color(colors[0])
    ax2 = ax.twinx()
    ss = pd.Series.from_csv(fsp, sep="\t")
    x = np.arange(ss.index[0], ss.index[-1], bin)
    ss = ss.values.reshape((len(ss) / bin, bin))
    ss = np.sum(ss, axis=1)
    ax2.plot(x, ss, color=colors[1], linewidth=1)
    ax2.set_ylabel("Subnucl. density")
    for t in ax2.get_yticklabels():
        t.set_color(colors[1])
    ax2.spines['right'].set_color(colors[1])
    ax2.spines['left'].set_color(colors[0])
    pylab.savefig(fout + ".pdf")


def plotProfile(fcn, fsp, fout):
    """
    MNase-seq profile plot, left y-axis is cn, right x-axis is sp
    """
    fig, ax = pylab.subplots(figsize=(2, 2.75 * 0.8))
    ax2 = ax.twinx()
    ss = pd.Series.from_csv(fcn, sep="\t")
    x = list(ss.index)
    y = list(smooth(np.array(ss.values)))
    x = np.arange(np.min(x), np.max(x), (np.max(x) - np.min(x)) / float(len(y)))
    ax.plot(x,y,color=colors[0],linewidth=1)
    #ax.plot(ss.index, ss, color=colors[0], linewidth=1)
    ax.set_ylabel("Nucleosome density")
    ax.set_xlabel("Distance from TSS")
    ax.set_title(fout)
    for t in ax.get_yticklabels():
        t.set_color(colors[0])
    ax2 = ax.twinx()
    ss = pd.Series.from_csv(fsp, sep="\t")
    x = list(ss.index)
    y = list(smooth(np.array(ss.values)))
    x = np.arange(np.min(x), np.max(x), (np.max(x) - np.min(x)) / float(len(y)))
    ax.plot(x,y,color=colors[1],linewidth=1)
    #ax2.plot(ss.index, ss, color=colors[1], linewidth=1)
    ax2.set_ylabel("Subnucl. density")
    for t in ax2.get_yticklabels():
        t.set_color(colors[1])
    ax2.spines['right'].set_color(colors[1])
    ax2.spines['left'].set_color(colors[0])
    fig.tight_layout()
    pylab.savefig(fout + ".pdf")


def main():
    gs = set(open("../../0.DHSs/ELIP_average.list").read().split("\n"))
    f = "../1.bedpe/WT_EILP.bedpe.gz"
    getCnSpProfiles(f,gs)
    gs = set(open("../../0.DHSs/ILCP_average.list").read().split("\n"))
    f = "../1.bedpe/WT_ILCP.bedpe.gz"
    getCnSpProfiles(f,gs)
    #fs = glob("../1.bedpe/*.bedpe.gz")
    #plotProfile( "data/KO_EILP_cn.txt","data/KO_EILP_sp.txt","KO_EILP" )
    #plotProfile( "data/WT_EILP_cn.txt","data/WT_EILP_sp.txt","WT_EILP" )
    #plotProfile( "data/WT_ILCP_cn.txt","data/WT_ILCP_sp.txt","WT_ILCP" )


if __name__ == '__main__':
    start_time = datetime.now()
    main()
    elapsed = datetime.now() - start_time
    print "The process is done"
    print "Time used:", elapsed
