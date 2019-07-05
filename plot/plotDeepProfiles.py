#--coding:utf-8--
"""
plotDeepProfiles.py
2019-06-03
"""

__author__ = "CAO Yaqiang"
__date__ = ""
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import os,gzip
from glob import glob

#3rd library
from tqdm import tqdm
import numpy as np
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
import pandas as pd
import brewer2mpl
colors = brewer2mpl.get_map('Set1', 'qualitative', 9).mpl_colors
colors.extend(brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors)
del colors[5]
del colors[10]



def getProfile(f,filterZero=True):
    ds = {}
    #-2000 and 2000 is the length before/after TSS
    xs = np.arange(-2000,2000,10)
    for i, line in enumerate(gzip.open(f)):
        if i == 0:
            continue
        line = line.split("\n")[0].split("\t")
        n = line[3]
        ss = list(map(float,line[6:]))
        if filterZero:
            if np.sum(ss) == 0.0:
            #if np.max(ss) < 0.2:
                continue
        #if line[5] == '-': do not need, as deeptools already consider it 
        #    ss.reverse()
        ss = pd.Series(ss,index=xs)
        ds[n] =ss
    ds = pd.DataFrame(ds).T
    return ds.mean()



def smooth(x, window_len=5, window='hanning'):
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


def profilePlot(fs,fout):
    fig, ax = pylab.subplots(figsize=(2,2.75*0.8))
    ax2 = ax.twinx()
    ss = getProfile(fs[0])
    x = list(ss.index)
    y = list(smooth(np.array(ss.values)))
    x = np.arange(np.min(x), np.max(x), (np.max(x) - np.min(x)) / float(len(y)))
    ax.plot(x,y,color=colors[0],linewidth=1)
    #ax.plot(ss.index,ss,color=colors[0],linewidth=1)
    ax.set_ylabel("Nucleosome density")
    ax.set_xlabel("Distance from TSS")
    ax.set_title(fout)
    for t in ax.get_yticklabels():
        t.set_color(colors[0])
    ax2 = ax.twinx()
    ss = getProfile(fs[1])
    x = list(ss.index)
    y = list(smooth(np.array(ss.values)))
    x = np.arange(np.min(x), np.max(x), (np.max(x) - np.min(x)) / float(len(y)))
    ax.plot(x,y,color=colors[1],linewidth=1)
    #ax2.plot(ss.index,ss,color=colors[1],linewidth=1)
    ax2.set_ylabel("Subnucl. density")
    for t in ax2.get_yticklabels():
        t.set_color(colors[1])
    ax2.spines['right'].set_color(colors[1])
    ax.spines['left'].set_color(colors[0])
    pylab.savefig(fout+".pdf") 


def main():
    profilePlot(["KO_EILP_cN_tss.txt.gz","KO_EILP_sP_tss.txt.gz"],"KO_EILP")
    profilePlot(["WT_EILP_cN_tss.txt.gz","WT_EILP_sP_tss.txt.gz"],"WT_EILP")
    profilePlot(["WT_ILCP_cN_tss.txt.gz","WT_ILCP_sP_tss.txt.gz"],"WT_ILCP")

main()
