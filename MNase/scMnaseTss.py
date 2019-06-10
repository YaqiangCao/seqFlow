#!/usr/bin/env python2.7
#--coding:utf-8--
"""
sepNcSpOther.py
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
from joblib import Parallel, delayed

#this
#from utils import getLogger, callSys, PET
from utils import *

#global settings
#logger
date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
logger = getLogger(fn=os.getcwd() + "/" + date.strip() + "_" +
                   os.path.basename(__file__) + ".log")



#TSS file
TSS = "/home/caoy7/caoy7/Projects/0.Reference/2.mm10/2.annotations/mm10_biomart_tss_tts_proteinCoding.bed"


def buildCovModel(readF,dfilter=[80,180]):
    """
    readF: bedpe.gz
    """
    print("building models for %s" % readF)
    n = readF.split('/')[-1].split(".bedpe.gz")[0]
    model = HTSeq.GenomicArray("auto", stranded=False)
    cn, sp, other = 0, 0, 0
    for i, line in enumerate(gzip.open(readF, 'rt')):
        line = line.split("\n")[0].split("\t")
        if len(line) < 6:
            continue
        if line[0] != line[3]:
            continue
        s = min(int(line[1]), int(line[4]))
        e = max(int(line[2]), int(line[5]))
        d = e - s
        if d >= dfilter[0] and d <= dfilter[1]: 
            m = (s+e)/2
            iv = HTSeq.GenomicInterval(line[0], m, m+1)
            model[iv] += 1
    return i, model


def getProfiles(f,ext=2000,bin=50,fo=None,skipZero=True):
    if fo == None:
        fo = "data/%s.txt"%(f.split("/")[-1].split(".bedpe")[0])
    if os.path.isfile(fo):
        return
    tot, model = buildCovModel(f)
    profile = np.zeros(2*ext/bin)
    for line in open(TSS):
        line = line.split("\n")[0].split("\t")
        m = int(line[1])
        s = m - ext 
        if s < 0:
            continue
        e = m + ext 
        iv = HTSeq.GenomicInterval(line[0],s,e,".")
        cvg = np.fromiter(model[iv], dtype='i', count=2*ext) 
        cvg = cvg.reshape((len(cvg)/bin,bin))
        #cvg = np.sum(cvg,axis=1)
        cvg = np.mean(cvg,axis=1)
        print(cvg)
        if skipZero and np.sum(cvg) == 0:
            continue
        if line[-1] == "-1":
            profile += cvg[::-1]
        else:
            profile += cvg
    s = pd.Series(profile,index=np.arange(0-ext,ext,bin))
    s.to_csv(fo,sep="\t")


def smooth(x, window_len=10, window='hanning'):
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


def plotSingleMnaseProfile(f,bin=20):
    ss = pd.Series.from_csv(f,sep="\t")
    x = list(ss.index)
    x = np.arange(x[0],x[-1],bin)
    y = ss.values
    y = y.reshape((len(y)/bin,bin))
    y = np.mean(y,axis=1)
    y = smooth(y,10)
    x = np.arange(np.min(x), np.max(x), (np.max(x) - np.min(x)) / float(len(y)))
    print(x.shape,y.shape)
    fig,ax = pylab.subplots()
    ax.plot(x,y)
    pylab.savefig(f.replace(".txt",".pdf"))


def plotMean():
    fs = glob("data/*.txt")
    ds = {}
    for f in fs:
        n = f.split("/")[-1].split(".")[0]
        s = pd.Series.from_csv(f,sep="\t")
        ds[n] = s
    ds = pd.DataFrame(ds)
    s = ds.mean(axis=1)
    fig,ax = pylab.subplots()
    ax.plot(s.index,s)
    pylab.savefig("test.pdf")


def main():
    #getProfiles("../../5.pool/1.bedpe/KO_EILP.bedpe.gz",bin=1)
    plotSingleMnaseProfile("data/KO_EILP.txt")
    #getProfiles("../../5.pool/1.bedpe/WT_EILP.bedpe.gz",bin=1)
    plotSingleMnaseProfile("data/WT_EILP.txt")

if __name__ == '__main__':
    start_time = datetime.now()
    main()
    elapsed = datetime.now() - start_time
    print "The process is done"
    print "Time used:", elapsed
