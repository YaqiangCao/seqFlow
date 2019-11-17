#!/usr/bin/env python2.7
#--coding:utf-8 --

"""
2015-12-29
2016-01-05: Modified as using highly used sets
"""


__author__="CAO Yaqiang"
__date__="2015-09-22"
__modified__=""
__email__="caoyaqiang0410@gmail.com"


#systematic library
import os,random
from glob import glob

#3rd 
import pandas as pd
import numpy as np

def callSys(cmds):
    """
    Call systematic commands without return.
    """
    for c in cmds:
         print c
         try:
             os.system(c)
         except:
             print "ERROR for %s" % c


def getfa(rs,ext,pre):
    bed = str(random.random())
    with open(bed,"w") as f:
        for i,nr in enumerate(rs):
            r = nr.split( "|" )
            if len(r) < 3:
                continue
            r[ 1 ] = int( r[ 1 ] )
            r[ 2 ] = int( r[ 2 ] )
            m = (r[1] + r[2])/2
            start = m - ext
            if start < 0 :
                start = 0
            end = m + ext
            line = [r[0],start,end,nr]
            f.write("\t".join(map(str,line))+"\n")
    cmd = "bedtools getfasta -fi /picb/molsysbio/usr/caoyaqiang/1.Projects/17.Alu/1.GenomeReference/1.hg38/1.hg38_Sequence/hg38.fa -bed {bed} -name -fo {fa}".format(bed=bed,fa=pre)
    callSys([cmd,"rm %s"%bed])



def getfas():
    for f in glob("../3.1.Modules/*.bed"):
        print os.path.split(f)[-1]
        ns = []
        rs = []
        for line in open(f):
            line = line.split( "\n" )[ 0 ].split( "\t" ) 
            s = int(line[2]) - int(line[1])
            ns.append(s)
            rs.append(line[-1])
        ns = np.array(ns)
        m,s = ns.mean(),int(ns.std())
        getfa(rs,3*s,os.path.split(f)[-1].replace(".bed",".fa"))


def main(  ):
    getfas()
    


if __name__=="__main__":
    main()
