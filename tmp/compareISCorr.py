#!/usr/bin/env python
#--coding:utf-8 --
"""
"""
__date__ = "2019-09-11"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#3rd library
import pandas as pd

#cLoops2
from cLoops2.settings import *


def gets(f):
    ds = {}
    for line in open(f):
        line = line.split("\n")[0].split("\t")
        ds[ "|".join( line[:-1]) ] = float( line[-1] )
    ds = pd.Series(ds)
    return ds

def main():
    sa = gets("trac1_IS.bdg")
    sb = gets("trac2_IS.bdg")
    sa = sa[sa<1]
    sb = sb[sb<1]
    ss = sa.index.intersection( sb.index)
    sa = sa[ss]
    sb = sb[ss]
    corr = sa.corr(sb)
    fig, ax = pylab.subplots()
    cmap = sns.cubehelix_palette(light=1, as_cmap=True)
    hb = ax.hexbin(sa,
                   sb,
                   gridsize=100,
                   cmap=cmap,
                   bins="log",
                   )
    cb = fig.colorbar(hb, ax=ax)
    cb.set_label('log10(N), number of points')
    ax.set_title("vector size:%s\nPCC:%.3f" % (len(sa), corr))
    ax.set_xlim([0.0,0.5])
    ax.set_ylim([0.0,0.5])
    ax.set_xlabel("Trac1 insulation score")
    ax.set_ylabel("Trac2 insulation score")
    pylab.savefig("Trac1Trac2_IS_scc.pdf")


if __name__ == "__main__":
    main()
