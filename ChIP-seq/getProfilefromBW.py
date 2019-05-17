#!/usr/bin/env python2.7
#--coding:utf-8--
"""
"""

__author__ = "CAO Yaqiang"
__date__ = "2015-04-06"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import glob, os, time
from datetime import datetime

#3rd library
import numpy as np
import pandas as pd
import bx.bbi.bigwig_file
from joblib import Parallel, delayed

#plotting settings
import matplotlib as mpl
mpl.use("pdf")
mpl.rcParams["pdf.fonttype"] = 42
import pylab, brewer2mpl
import seaborn as sns
colors = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors
colors.extend(brewer2mpl.get_map('Set3', 'qualitative', 8).mpl_colors)
sns.set(style="white")


def smooth(x, window_len=11, window='hanning'):
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
    #print(len(s))
    if window == 'flat':  #moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.' + window + '(window_len)')

    y = np.convolve(w / w.sum(), s, mode='valid')
    return y


def getProfile(bed, bw, pre=None, ext=5000, shift=0, bins=200):
    if pre == None:
        pre = bw.split("/")[-1].split(".")[0]
    data = {}
    bwh = bx.bbi.bigwig_file.BigWigFile(open(bw, "rb"))
    x = np.arange(ext * 2) - ext
    x = x.reshape(
        -1,
        bins,
    ).mean(axis=1)
    for line in open(bed):
        line = line.split("\n")[0].split("\t")
        submit = int(line[1]) + shift
        start = submit - ext
        if start < 0:
            continue
        end = submit + ext
        d = bwh.get_as_array(line[0], start, end)
        for i in xrange(len(d)):
            if np.isnan(d[i]):
                d[i] = 0.0
        #get bins
        nd = d.reshape(
            -1,
            bins,
        ).mean(axis=1)
        data[line[-1]] = nd
    data = pd.DataFrame(data).T
    data.columns = x
    data.to_csv(pre + ".txt", sep="\t")


def plotProfile(pre):
    fig, ax = pylab.subplots()
    fs = glob.glob("*.txt")
    fs.sort()
    for i, f in enumerate(fs):
        label = f.split(".")[0]
        mat = pd.read_table(f, index_col=0)
        x = mat.columns
        x = map(float, x)
        y = mat.mean(axis=0)
        y = smooth(y)
        x = np.arange(np.min(x), np.max(x),
                      (np.max(x) - np.min(x)) / float(len(y)))
        ax.plot(x, y, label=label, color=colors[i], lw=3)
    leg = ax.legend(loc="upper left", fancybox=True, bbox_to_anchor=(1, 1))
    ax.set_xlabel("Distance to TSS (bp)")
    ax.set_ylabel("Smoothed Normalized Reads Density")
    pylab.savefig(pre + ".pdf", bbox_inches="tight")


def plotProfile2(pre):
    """
    Minus input signal.
    """
    fig, ax = pylab.subplots()
    fs = glob.glob("*.txt")
    fs.sort()
    ds = {}
    for i, f in enumerate(fs):
        label = f.split(".")[0]
        mat = pd.read_table(f, index_col=0)
        x = mat.columns
        x = map(float, x)
        y = mat.mean(axis=0)
        s, t = label.split("_")
        if s not in ds:
            ds[s] = {"T": 0, "I": 0, "S": t}
        if t == "input":
            ds[s]["I"] = [x, y]
        else:
            ds[s]["T"] = [x, y]
    for i, key in enumerate(ds.keys()):
        x = ds[key]["T"][0]
        y = ds[key]["T"][1] - ds[key]["I"][1]
        y = smooth(y)
        x = np.arange(np.min(x), np.max(x),
                      (np.max(x) - np.min(x)) / float(len(y)))
        label = key + "_" + ds[key]["S"]
        ax.plot(x, y, label=label, color=colors[i], lw=3)
    leg = ax.legend(loc="upper left", fancybox=True, bbox_to_anchor=(1, 1))
    ax.set_xlabel("Distance to TSS (bp)")
    ax.set_ylabel("Smoothed Normalized Reads Density")
    pylab.savefig(pre + ".pdf", bbox_inches="tight")


def main():
    bed = "../AW_Enhanced.bed"
    bws = "/picb/molsysbio/usr/caoyaqiang/1.Projects/2.WL_hESC/2.ChIPSeq/18.Selh/1.SMAD2/4.Pileups"
    bws = glob.glob(bws + "/*.bw")
    Parallel(n_jobs=len(bws))(delayed(getProfile)(bed, bw) for bw in bws)
    plotProfile2("SMAD2")
    #plotProfile( "SMAD2" )


main()
