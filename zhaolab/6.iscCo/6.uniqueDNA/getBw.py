
import os
import time
from glob import glob

#3rd
import click
from tqdm import tqdm
import pyBigWig
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

#pipi
from pipi.sets import *
from pipi.utils import getLogger, callSys, isTool

#global settings
#logger
date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
logger = getLogger(fn=os.getcwd() + "/" + date.strip() + "_" +
                   os.path.basename(__file__) + ".log")



def getTss(bw,tssf,ext=2500,skipZeros=True):
    n = bw.split("/")[-1].split(".bw")[0]
    print(n,tssf)
    bw = pyBigWig.open(bw)
    x = range(-ext,ext+1)
    c = 0
    ys = np.zeros(len(x))
    for line in tqdm(open(tssf).read().split("\n")):
        line = line.split("\n")[0].split("\t")
        if len(line) < 3:
            continue
        chrom = line[0]
        center =int( (int(line[1]) + int(line[2] ))/2 ) 
        start = center -ext 
        end = center + ext+1
        if start < 0:
            continue
        try:
            s = np.array(bw.values(chrom, start, end))
        except:
            continue
        s = np.nan_to_num(s)
        if line[5] == "-":
            s = list(s)
            s.reverse()
            s = np.array(s)
        if skipZeros:
            if np.sum(s) > 0:
                ys = ys + s
                c += 1
        else:
            ys = ys + s
            c += 1
    ys = ys / c
    fig, ax = pylab.subplots()
    ax.plot(x, ys)
    ax.set_xlabel("Distance to TSS (bp)") 
    ax.set_ylabel("RPM")
    pylab.savefig(n+"_tss.pdf")


@click.command()
@click.option("-name", required=True, help="Sample name/id for the data.")
@click.option("-org",
              required=True,
              help="Organism for the data.",
              type=click.Choice(["hg38", "mm10"]))
def main(name,org):
    genomeRef = {
        "hg38":
        "/home/caoy7/caoy7/Projects/0.Reference/1.hg38/1.fa/hg38.chrom.sizes",
        "mm10":
        "/home/caoy7/caoy7/Projects/0.Reference/2.mm10/1.fa/mm10.chrom.sizes",
    }
    tssRef = {
        "hg38":"/home/caoy7/caoy7/Projects/0.Reference/1.hg38/2.annotations/gencode_v30_pcRNA_tss.bed",
        "mm10":"/home/caoy7/caoy7/Projects/0.Reference/2.mm10/2.annotations/gencode_vM21_pcRNA_tss.bed",
        }
    bed = name + "_DNA.bed"
    bam = name +"_DNA.bam"
    bw = name+"_DNA.bw"
    c1 = "bedtools bedtobam -i {bed} -g {genome} > {bam}".format( bed=bed, genome=genomeRef[org], bam=bam )
    c2 = "samtools sort -@ 2 -o {bam} {bam}".format(bam=bam)
    c3 = "samtools index {bam} {bai}".format(bam=bam,bai=bam.replace(".bam",".bai"))
    c4 = "bamCoverage -b {bam} -o {bw} -p 10 --normalizeUsing CPM".format(bam=bam,bw=bw)
    callSys([c1,c2,c3,c4], logger)
    print(tssRef[org])
    getTss( bw, tssRef[org] )


main()
