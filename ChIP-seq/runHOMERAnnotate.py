#!/usr/bin/env python3.6
#--coding:utf-8 --
"""
runHomerAnnotate.py
2019-05-13: updated
"""

#sys
import os, subprocess
from glob import glob
from datetime import datetime

#3rd
import pandas as pd
from joblib import Parallel, delayed


def runHOMER(bed, genome="mm10"):
    name = bed.split("/")[-1].split(".")[0]
    stat = name + "_stat.txt"
    ano = name + "_ano.txt"
    if os.path.exists(ano):
        print(ano, "generated, return.")
        return
    cmd = "annotatePeaks.pl {bed} {genome} -go {go_dir} -annStats {stat} > {ano}".format(
        bed=bed,
        genome=genome,
        stat=stat,
        go_dir="go_" + name,
        ano=ano)
    try:
        print(cmd)
        subprocess.call(cmd, shell=True)
    except:
        return


def parseTarget():
    fs = glob("*ano.txt")
    for f in fs:
        print(f)
        pre = f.split(".")[0]
        mat = pd.read_table(f, index_col=0)
        gs = {}
        for i in mat.index:
            g = mat.loc[i, "Nearest PromoterID"]
            if g not in gs:
                gs[g] = set()
            gs[g].add(i)
        with open(pre + ".gmt", "w") as nf:
            for g, rs in gs.items():
                line = g + "\t" + ",".join(rs) + "\n"
                nf.write(line)


def main():
    beds = glob("../../1.peaks/*.narrowPeak")
    Parallel( n_jobs=1 )( delayed( runHOMER )( f, "mm10") for f in beds )
    #parseTarget()


if __name__ == "__main__":
    startTime = datetime.now()
    main()
    elapsed = datetime.now() - startTime
    print("The process is done.")
    print("Time used:", elapsed)
