#!/usr/bin/env python2.7
#--coding:utf-8--
"""
run_HOMER_profile.py
2015-03-26: Specific plot as minus input.
"""

__author__ = "CAO Yaqiang"
__date__ = "2015-03-05"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import glob, os, time
from datetime import datetime

#3rd library
import matplotlib as mpl
mpl.use("pdf")
mpl.rcParams["pdf.fonttype"] = 42
import pylab, brewer2mpl
import numpy as np
import seaborn as sns
import pandas as pd
import brewer2mpl
colors = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors
colors.extend(brewer2mpl.get_map('Set3', 'qualitative', 8).mpl_colors)
sns.set(style="white")

#own library


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


def getProfile(gs, term, org="hg19", tags=None):
    tags = glob.glob(tags + "/*%s*" % term)
    for tag in tags:
        if not os.path.isdir(tag):
            continue
        pre = tag.split("/")[-1]
        cmd = "annotatePeaks.pl tss {org} -list {gs} -size 10000 -hist 50 -ghist -d {tag} > {pre}.txt".format(
            org=org, gs=gs, tag=tag, pre=pre)
        print cmd
        os.system(cmd)


def plotProfile(pre):
    fig, ax = pylab.subplots()
    fs = glob.glob("*.txt")
    fs.sort()
    for i, f in enumerate(fs):
        label = f.split(".")[0]
        mat = pd.read_table(f, index_col=0)
        x = mat.columns
        x = map(int, x)
        y = mat.mean(axis=0)
        y = smooth(y)
        x = np.arange(
            np.min(x), np.max(x), (np.max(x) - np.min(x)) / float(len(y)))
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
        x = map(int, x)
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
        x = np.arange(
            np.min(x), np.max(x), (np.max(x) - np.min(x)) / float(len(y)))
        label = key + "_" + ds[key]["S"]
        ax.plot(x, y, label=label, color=colors[i], lw=3)
    #leg = ax.legend( loc="upper left",fancybox=True,bbox_to_anchor=( 1,1 ) )
    leg = ax.legend(loc="best", fancybox=True)
    ax.set_xlabel("Distance to TSS (bp)")
    ax.set_ylabel("Smoothed Normalized Reads Density")
    pylab.savefig(pre + ".pdf", bbox_inches="tight")


def main():
    gs = "../../AW_Enhanced.list"
    term = ""
    tags = "../../../../../2.ChIPSeq/18.Selh/1.SMAD2/3.Tags/"
    #getProfile( gs,term,"hg19",tags )
    plotProfile2("SMAD2")
    #plotProfile( "SMAD2" )


main()
