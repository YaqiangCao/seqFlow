
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
    gtfRef = {
        "hg38":"/home/caoy7/caoy7/Projects/0.Reference/1.hg38/2.annotations/gencode_human_v30_exons.gtf",
        "mm10":"/home/caoy7/caoy7/Projects/0.Reference/2.mm10/2.annotations/gencode_v21_exons.gtf",
        }
    bam = name +"_RNA.bam"
    bw = name+"_RNA.bw"
    gtf = gtfRef[org]
    c1 = "samtools sort -@ 2 -o {bam} {bam}".format(bam=bam)
    c2 = "samtools index {bam} {bai}".format(bam=bam,bai=bam.replace(".bam",".bai"))
    c3 = "bamCoverage -b {bam} -o {bw} -p 10 --normalizeUsing CPM".format(bam=bam,bw=bw)
    c4 = "/home/caoy7/caoy7/Packages/qualimap_v2.2.1/qualimap rnaseq -bam {bam} -outdir qc -a proportional -gtf {gtf}".format(bam=bam,gtf=gtf)
    callSys([c1,c2,c3,c4], logger)


main()
