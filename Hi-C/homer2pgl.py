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



def homer2pgl(fin,fout):
    if os.path.isfile(fout):
        return
    cs = ["chr%s"%i for i in xrange(1,23)]
    cs.extend(["chrM","chrX","chrY"])
    with open(fout,"w") as fo:
        for i,line in enumerate(open(fin)):
            if i == 0:
                continue
            line = line.split( "\n" )[ 0 ].split( "\t" ) 
            if line[2] != line[8]:
                continue
            if line[2] not in cs:
                continue
            p = float(line[17])
            fdr = float(line[18])
            #print line, p, fdr
            if fdr < 0.05:
                if int(line[3]) < 0 or int(line[9]) <0:
                    continue
                nline = [line[2],line[3],line[4],line[8],line[9],line[10],line[0],line[17],line[14]]
                fo.write("\t".join(map(str,nline))+"\n")


for f in glob("./*.txt"):
    homer2pgl(f,f.split("/")[-1].replace(".txt","_HOMER.pgl"))

