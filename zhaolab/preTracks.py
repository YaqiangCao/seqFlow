#!/usr/bin/env python
#--coding:utf-8--
"""
preTracks.py
2019-06-28: Prepare tracks to upload to washU.
"""

#sys
import os, time
from glob import glob
from datetime import datetime

#3rd
import pandas as pd

#seqFlow
from utils import isTool, getLogger, callSys

#global setting
#logger
date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
logger = getLogger(fn=os.getcwd() + "/" + date.strip() + "_" +
                   os.path.basename(__file__) + ".log")





def main():
    metaf = "../1.fastq/1_20190628_KZ1848_KZ1841_fq.txt"
    root = "/var/www/html/browser/data/mm10/CYQ/"
    meta = pd.read_csv(metaf,sep="\t",index_col=0)
    ds = {}
    with open("records.txt","w") as fo:
        line = ["BEDgraph","type","Scientist","Cell-type","Assay-type","Description"]
        fo.write("\t".join(line)+"\n")
        for g in meta.index:
            t = g + "_washU.txt.gz"
            if not os.path.isfile(t):
                print("%s not exist for %s"%(t,g))
                continue
            ng = root + t
            ct = "unkwn"
            nline = [ng,"bedGraph","CYQ-UNKWN",ct,"CHIC",g+"_"+meta.loc[g,"Sample discription"]]
            fo.write("\t".join(nline)+"\n")


if __name__ == "__main__":
    startTime = datetime.now()
    main()
    elapsed = datetime.now() - startTime
    print("The process is done.")
    print("Time used:", elapsed)
