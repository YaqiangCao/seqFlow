#sys
import os, random, gzip
from glob import glob
from collections import Counter

#3rd library
import HTSeq
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
import seaborn as sns
import brewer2mpl
colors =  brewer2mpl.get_map( 'Set3', 'qualitative',8 ).mpl_colors
pgltool="/picb/molsysbio/usr/caoyaqiang/4.ENV/SG/sage-7.3-cyq/upstream/bio/pgltools/sh/pgltools"



def sortPgl(f):
    if os.path.isfile("%s.sorted"%f):
        return
    cmd = "%s sort %s > %s.sorted"%(pgltool,f,f)
    print cmd
    os.system(cmd)


def pgl2washU(fin):
    fout = fin.replace(".pgl.sorted",".washu.txt")
    with open(fout,"w") as fo:
        for line in open(fin):
            line = line.split( "\n" )[ 0 ].split( "\t" ) 
            t1 = "%s:%s-%s"%(line[0],line[1],line[2])
            t2 = "%s:%s-%s"%(line[3],line[4],line[5])
            nline = [t1,t2,"1"]
            fo.write("\t".join(nline)+"\n")


def pgl2juice(fin):
    fout = fin.replace(".pgl.sorted",".juicebox.txt")
    if os.path.isfile(fout):
        return
    cmd = "%s juicebox -a %s > %s"%(pgltool,fin,fout)
    print cmd
    os.system(cmd)

"""
for f in glob("*.pgl"):
    sortPgl(f)
"""

for f in glob("*.sorted"):
    pgl2washU(f)
    #pgl2juice(f)

