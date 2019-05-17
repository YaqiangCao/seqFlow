#!/usr/bin/env python2.7
#--coding:utf-8--
"""
Scan all known motif appearance nearby peaks, than plot
2015-03-26: Modified result parse and plot
"""

__author__ = "CAO Yaqiang"
__date__ = "2014-05-25"
__modified__ = "2015-03-26"
__license__ = "GPL"
__version__ = "0.1"
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import glob, os
from datetime import datetime

#3rd library
import matplotlib as mpl
mpl.use("pdf")
mpl.rcParams["pdf.fonttype"] = 42
import pylab, brewer2mpl
import numpy as np
#import seaborn as sns
import pandas as pd
import brewer2mpl
colors = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors
colors.extend(brewer2mpl.get_map('Set3', 'qualitative', 8).mpl_colors)

#sns.set( style="white" )


def homerMotifScan(bed, genome="hg19"):
    fn = bed.split("/")[-1].split(".")[0] + ".txt"
    cmd = "annotatePeaks.pl {bed} {genome} -hist 20 -size 10000 -m /home/caoyaqiang/caoyaqiang_han/4.ENV/2.Bios/HOMER/data/knownTFs/vertebrates/known.motifs > {fn}".format(
        bed=bed, genome=genome, fn=fn)
    print cmd
    os.system(cmd)


def parseMotifs(f):
    mat = pd.read_table(f, index_col=0)
    mat = mat.iloc[:, :-4]
    mat = mat.iloc[:, 0::3]
    cs = []
    for c in mat.columns:
        c = c.split(" ")[0].split("/")
        c[0] = c[0].split("(")[0]
        if len(c) > 1 and "GSE" in c[1]:
            c[1] = c[1].split("(")[1][:-1]
        c = "_".join(c)
        cs.append(c)
    mat.columns = cs
    fn = f.split("/")[-1].split(".")[0] + "_parsed.txt"
    mat.to_csv(fn, sep="\t")


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


def plotMotifs(
        f,
        ts=[
            "ZBTB33_GSE32465_Homer",
            "Smad2_GSE29422_Homer",
        ],
):
    pre = f.split("_")[0]
    mat = pd.read_table(f, index_col=0)
    x = map(int, mat.index)
    fig, ax = pylab.subplots()
    for i, t in enumerate(ts):
        y = mat[t].values
        y = smooth(y)
        x = np.arange(np.min(x), np.max(x),
                      (np.max(x) - np.min(x)) / (float(len(y))))
        label = t.split("_")[0]
        ax.plot(x, y, label=label, color=colors[i], lw=2)
    #leg = ax.legend( loc="upper left",fancybox=True,bbox_to_anchor=( 1,1 ) )
    leg = ax.legend(loc="best", fancybox=True)
    ax.set_xlabel("Distance to Peaks Summits (bp)")
    ax.set_ylabel("Motifs Density")
    ax.set_title(pre)
    pylab.savefig(pre + ".pdf", bbox_inches="tight")


#beds = glob.glob( "../1.summits/*SMAD2*" )
#beds.extend( glob.glob( "../1.summits/*b-catenin*" ) )
#map( homer_known_motif_scan,beds )
#fs = glob.glob("*.txt")
#map( parseMotifs,fs )
fs = glob.glob("*parsed.txt")
map(plotMotifs, fs)
