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





def loops2degreesSharp(fin,fout):
    model = HTSeq.GenomicArrayOfSets("auto",stranded=0)
    for i, line in enumerate(open(fin)):
        if i == 0:
            continue
        line = line.split( "\n" )[ 0 ].split( "\t" )  
        iva = HTSeq.GenomicInterval(line[0],int(line[1]),int(line[2]))
        ivb = HTSeq.GenomicInterval(line[3],int(line[4]),int(line[5]))
        model[iva] += line[6]+"-left"
        model[ivb] += line[6]+"-right"
    with open(fout,"w") as fo:
        for iv, value in list(model.steps()):
            if value == set([]):
                continue
            line = [iv.chrom,iv.start,iv.end,len(value),",".join(value)]
            fo.write("\t".join(map(str,line))+"\n")


def loops2degreesBroad(fin,fout):
    if os.path.isfile(fout):
        return
    model = HTSeq.GenomicArrayOfSets("auto",stranded=0)
    model2 = HTSeq.GenomicArrayOfSets("auto",stranded=0)
    for i, line in enumerate(open(fin)):
        line = line.split( "\n" )[ 0 ].split( "\t" )  
        iva = HTSeq.GenomicInterval(line[0],int(line[1]),int(line[2]))
        ivb = HTSeq.GenomicInterval(line[3],int(line[4]),int(line[5]))
        try:
            model[iva] += 1
            model[ivb] += 1 
            model2[iva] += line[6]+"-left"
            model2[ivb] += line[6]+"-right"
        except:
            print fin,line
    with open(fout,"w") as fo:
        for iv, value in list(model.steps()):
            if value == set([]):
                continue
            ds = set()
            for ivb,valueb in model2[iv].steps():
                ds.update(valueb)
            ds = list(ds)
            line = [iv.chrom,iv.start,iv.end,len(ds),",".join(ds)]
            fo.write("\t".join(map(str,line))+"\n")



for f in glob("../1.sets/*.pgl.sorted"):
    #loops2degreesSharp(f,f.split("/")[-1].replace(".pgl.sorted",".sharpDegree"))
    loops2degreesBroad(f,f.split("/")[-1].replace(".pgl.sorted",".broadDegree"))

