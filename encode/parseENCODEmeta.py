import os,sys
from glob import glob
import pandas as pd


def parseMeta(fs,metaf) :
    """
    Rename ENCODE tracks
    """
    meta = pd.read_csv(metaf,sep="\t",index_col=0)
    for f in fs:
        n = f.split("/")[-1].split(".")[0]
        ft = f.split("/")[-1].split(".")[-1]
        if n not in meta.index:
            continue
        c = meta.loc[n,"Biosample term name"]
        t = meta.loc[n,"Experiment target"].split("-")[0]
        br = meta.loc[n,"Biological replicate(s)"]
        tr = meta.loc[n,"Technical replicate"]
        nn = "%s_%s_bio%s"%(c,t,br)
        if not pd.isnull(tr):
            nn += "_tech%s"%tr
        nn = nn+"."+ft
        cmd = "ln -s %s %s"%(f,nn)
        print(cmd)
        os.system(cmd)

fs = glob("../1.bam/*.bam")
metaf = "../1.bam/metadata.tsv"
parseMeta( fs, metaf)
