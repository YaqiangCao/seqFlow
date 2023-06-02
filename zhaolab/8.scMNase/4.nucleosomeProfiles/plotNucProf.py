#!/usr/bin/env python
#--coding:utf-8--
"""
plotNucProf.py 
Plot nucleosome and subnucleosome profiles around TF binding/DHSs summits.
"""

__author__ = "CAO Yaqiang"
__date__ = "2023-04-06"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#sys
import os
import gzip
import subprocess
from datetime import datetime

#3rd library
import click
import pyBigWig
import numpy as np
import pandas as pd
from tqdm import tqdm
#mainly used for plotting settings
from cLoops2.settings import *


def getSignal(bw, rs, ext=1000, skipZeros=True):
    x = range(-ext, ext + 1)
    c = 0
    ys = np.zeros(len(x))
    bw = pyBigWig.open(bw)
    for r in tqdm(rs):
        try:
            s = np.array(bw.values(r[0], r[1], r[2]))
        except:
            continue
        s = np.nan_to_num(s)
        if skipZeros:
            if np.sum(s) > 0:
                ys = ys + s
                c += 1
        else:
            ys = ys + s
            c += 1
    ys = ys / c
    return ys


def plotNucProfile(fc, fn, fs, fo, ext=1000, title=""):
    """
    Plot nucleosome and subnucleosome profiles around some centers.
    """
    #read all the centers
    rs = []
    for line in open(fc).read().split("\n"):
        line = line.split("\n")[0].split("\t")
        if len(line) < 2:
            continue
        chrom = line[0]
        center = int(line[1])
        start = center - ext
        end = center + ext + 1
        if start < 0:
            continue
        rs.append([chrom, start, end])
    rs = sorted(rs, key=lambda x: x[0])

    #get signals around centers
    x = range(-ext, ext + 1)
    report = "\t".join([str(datetime.now()), f"Getting signals from {fn}"])
    nys = getSignal(fn, rs, ext)
    report = "\t".join([str(datetime.now()), f"Getting signals from {fs}"])
    sys = getSignal(fs, rs, ext)

    #plot signals
    fig, ax = pylab.subplots(figsize=(3.2, 2.2))
    ax.plot(x, nys, color=colors[0], linewidth=1)
    #ax.set_ylabel("Nucleosome center-weighted occupancy")
    ax.set_ylabel("Nucleosome density")
    ax.yaxis.label.set_color(colors[0])
    ax.set_xlabel("Distance from center (bp)")
    ax.spines['left'].set_color(colors[0])
    for t in ax.get_yticklabels():
        t.set_color(colors[0])
    ax2 = ax.twinx()
    ax2.plot(x, sys, color=colors[1], linewidth=1)
    ax2.set_ylabel("Subnucleosome density")
    ax2.yaxis.label.set_color(colors[1])
    for t in ax2.get_yticklabels():
        t.set_color(colors[1])
    ax2.spines['right'].set_color(colors[1])
    ax2.spines['left'].set_color(colors[0])
    ax.set_title(title)
    #sns.despine(ax=ax, top=True, right=False, left=False, bottom=False, offset=0, trim=True)
    #sns.despine(ax=ax2, top=True, right=False, left=False, bottom=False, offset=0, trim=True)
    #fig.tight_layout()
    pylab.savefig(fo + "_NucProf.pdf")


@click.command()
@click.option(
    "-c",
    required=True,
    help=
    "Input center file in .txt format. First column is chromosome and second column is coordinate,seperated by \\t",
    type=str,
)
@click.option(
    "-n",
    required=True,
    help=
    "Input nucleosome signals in bigWig format, as _cwos.bw file outputed by getNucProf.py.",
    type=str,
)
@click.option(
    "-s",
    required=True,
    help=
    "Input subnucleosome signals (or other signal as TF bindings) in bigWig format, as _sub.bw file outputed by getNucProf.py.",
    type=str,
)
@click.option(
    "-o",
    required=True,
    help="Output file prefix",
    type=str,
)
@click.option(
    "-ext",
    required=False,
    default=1000,
    help="Extension from the center. Default is 1000 bp.",
    type=int,
)
@click.option(
    "-title",
    required=False,
    default="",
    help=
    "Title for the figure to indicate the center type, such as DHSs, CTCF summits.",
    type=str,
)
def run(c, n, s, o, ext=1000, title=""):
    """
    Plot nucleosome and subnucleosome profiles around centers (non-directional) using the signals output by getNucProf.py

    Example: 

    plotNucProf.py -c naive_ctcf_summits.txt -n naive_cwos.bw -s naive_subn.bw -o naive
    """
    #start
    start = datetime.now()
    cmd = os.path.basename(__file__)
    report = "\t".join([
        str(datetime.now()),
        f"{cmd} -c {c} -n {n} -s {s} -o {o} -ext {ext} -title {title}"
    ])
    print(report)

    plotNucProfile(c, n, s, o, ext=ext, title=title)

    #finished
    end = datetime.now()
    usedTime = end - start
    report = "\t".join(
        [str(datetime.now()), f"{cmd} job finished. Used time: {usedTime}"])
    print(report + "\n\n")


if __name__ == "__main__":
    run()
