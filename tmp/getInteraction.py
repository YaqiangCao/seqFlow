#!/usr/bin/env python2.7
#--coding:utf-8--

#sys
import os
from glob import glob

#3rd library
import HTSeq
import pandas as pd


def readRep( f ):
    repModel = HTSeq.GenomicArrayOfSets( "auto", stranded = False )
    for line in open( f ):
        line = line.split( "\n" )[ 0 ].split( "\t" ) 
        locus = [ line[ 0 ],line[ 1 ],line[ 2 ],line[ 3 ] ]
        try:
            iv = HTSeq.GenomicInterval( locus[ 1 ],int( locus[ 2 ] ),int( locus[ 3 ] ) )
            repModel[ iv ] += locus[0]
        except:
            continue
    return repModel


def getRepOverlap(repModel,iv):
    reps = set(  )
    for niv, value in list( repModel[ iv ].steps(  ) ):
        reps.update( value  )
    return reps


def getAllInt():
    repModel = readRep("/picb/molsysbio/usr/caoyaqiang/1.Projects/17.Alu/1.GenomeReference/1.hg38/6.Repeats_Annotation/1.GenomicLocation/hg38Rep.txt")
    for fin in glob("../1.loops/*.pgl"):
        print fin
        fout = fin.split("/")[-1]
        if os.path.isfile(fout):
            continue
        with open(fout,"w") as f:
           for i,line in enumerate(open(fin)):
               if i == 0:
                   continue
               line = line.split( "\n" )[ 0 ].split( "\t" ) 
               a = HTSeq.GenomicInterval(line[0],int(line[1]),int(line[2])) 
               b = HTSeq.GenomicInterval(line[3],int(line[4]),int(line[5])) 
               ar = getRepOverlap(repModel,a)
               br = getRepOverlap(repModel,b)
               if len(ar) == 0 and len(br) == 0:
                   continue
               line.extend([",".join(ar),",".join(br)])
               f.write("\t".join(line)+"\n")


if __name__ == '__main__':
    getAllInt()
    
