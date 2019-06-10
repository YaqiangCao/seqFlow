#!/usr/bin/env python2.7
#--coding:utf-8--


"""
addFeature.py
"""



__author__="CAO Yaqiang"
__date__=""
__modified__=""
__email__="caoyaqiang0410@gmail.com"


#systematic library

#3rd
import pandas as pd


def addFeature( fa,fb,c ):
    mat = pd.read_table( fa,index_col=0 )
    old2new = {}
    for i in mat.index:
        old2new["|".join(i.split("|")[:-1])] = i
    old2new = pd.Series(old2new)
    s = pd.Series.from_csv( fb,sep="\t" )
    s.index = old2new[s.index]
    s = s[ mat.index ]
    mat[ c ] = s
    mat.to_csv( fa,sep="\t" )


def main(  ):
    fa = "mm10Rep.txt"
    #add mappability
    fb = "../2.mappability/mappabilityS36.txt"
    c = "mappability_36"
    addFeature( fa,fb,c )
   # fb = "../2.mappability/mappabilityS50.txt"
   # c = "mappability_50"
   # addFeature( fa,fb,c )
    fb = "../3.GC/mm10Rep.txt"
    c = "GC"
    addFeature( fa,fb,c )


main(  )
