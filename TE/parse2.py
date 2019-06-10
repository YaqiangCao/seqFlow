#!/usr/bin/env python2.7
#--coding:utf-8--


"""
parse2.py
"""



__author__="CAO Yaqiang"
__date__=""
__modified__=""
__email__="caoyaqiang0410@gmail.com"


#systematic library

#3rd
import pandas as pd


def readRepFamily(  ):
    f = "/picb/molsysbio/usr/caoyaqiang/1.Projects/17.Alu/1.GenomeReference/3.mm10/3.RepeatsAnnotation/repFamily.gmt"
    #read family info
    rep2Family = {  }
    for line in open( f ):
        line = line.split( "\n" )[ 0 ].split( "\t" )  
        reps = line[ 1 ].split( "," )
        for rep in reps:
            rep2Family[ rep ] = line[ 0 ]
    return rep2Family


def getBg( ext=0 ):
    rep2fam = readRepFamily(  )
    f = "mm10_rep_locus.txt"
    ds = {  }
    with open(  "mm10Rep.txt","w" ) as nf:
        nf.write( "rep\tchr\ts\te\tstrand\tsubfamily\tfamily\tlocus\n" )
        for i,line in enumerate( open( f ) ):
            if i % 10000 == 0:
                print i
            line = line.split( "\n" )[ 0 ].split( "\t" ) 
            rep = line[ 4 ]
            if rep not in rep2fam:
                continue
            repf = rep2fam[ rep ]
            rep = "|".join( line )
            if ext == 0:
                s = int( line[ 1 ] )
                e = int( line[ 2 ] )
            else:
                c = 1/2.0*( int( line[ 1 ] )+int( line[ 2 ] ) )
                s = c - ext
                if s < 0:
                    s = 0
                e = c +  ext
            nline = [ rep,line[ 0 ],line[ 1 ],line[ 2 ],line[ 3 ],line[ 4 ],repf,line[ 5 ]]
            nf.write( "\t".join( map( str,nline ) ) + "\n" )

getBg(  )
