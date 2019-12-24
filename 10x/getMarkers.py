import os
from glob import glob

#3rd
import numpy as np
import pandas as pd

def parseMarkers(f):
    """
    Prase the markers in GMT format. 
    """
    markers = {}
    for line in open(f):
        line = line.split("\n")[0].split("\t")
        markers[ line[0] ] = set(  line[1].split(",") )
    return markers

def getMarkerExp(mat,gs):
    ns = []
    gs = set([g.upper() for g in list(gs)])
    for i in mat.index:
        n = i.split("|")[1].upper()
        if n in gs:
            ns.append( i )
    return mat.loc[ns,]
        

def get(f,markers):
    mat = pd.read_csv(f,index_col=0,sep=",")
    ms = []
    for i,k in enumerate(markers.keys()):
        nk = k.replace(" ","_")
        m = getMarkerExp( mat, markers[k])
        m.index = [ str(i)+"|"+j+"|"+nk for j in m.index]
        ms.append( m )
    mat = pd.concat(ms)
    n = f.split("/")[-1]
    mat.to_csv(n.replace(".txt","_marker.txt"),sep="\t")

def main():
    f = "../3.matrix/GR_1.txt"
    markers = parseMarkers("../0.reference/mouse_spleen_jian_marker.txt")
    for f in glob("../3.matrix/*.txt"):
        get(f,markers)
    """
    for ka in markers.keys():
        for kb in list(markers.keys()):
            if ka == kb:
                continue
            print(ka,kb, markers[ka].intersection( markers[kb]))
    """

main()
