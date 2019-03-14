import os,commands,shutil,gzip
from datetime import datetime
from glob import glob
import pandas as pd

def runsys(cmd):
    print cmd
    s = datetime.now()
    #_code,output = commands.getstatusoutput(cmd)
    os.system(cmd)
    e = datetime.now()
    t = e-s
    return t


def runFitHiC():
    chrs = [ "chr%s"%i for i in xrange(1,22)]
    chrs.extend(["chrM","chrX","chrY"])
    totalTime = 0
    for i,c in enumerate(chrs):
        fbias = "../sep/%s_fithic.biases.gz"%c
        fmap = "../sep/%s_fithic.fragmentMappability.gz"%c
        fi = "../sep/%s_fithic.interactionCounts.gz"%c
        cmd = "fithic -f {fmap} -i {fi} -t {fbias} -o ./ -q -l {c} -L 20000 ".format(fmap=fmap,fi=fi,fbias=fbias,c=c)
        t = runsys(cmd)
        if i == 0:
            totalTime = t
        else:
            totalTime = totalTime + t
    print totalTime


def summaryFitHic(qcut=0.05):
    fs = glob("*.gz")
    fs.sort()
    ds = set()
    with open("fithic_summary.txt","w") as fo:
        for f in fs:
            print f
            for i,line in enumerate(gzip.open(f)):
                line = line.split( "\n" )[ 0 ].split( "\t" ) 
                if i == 0:
                    continue
                if len(line) < 7:
                    continue
                q = float(line[-1])
                if q > qcut:
                    continue
                key = (line[0],line[1],line[3])
                if key not in ds:
                    fo.write("\t".join(line)+"\n")
                    ds.add(key)



runFitHiC()
summaryFitHic()
