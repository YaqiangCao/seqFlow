#!/usr/bin/env python
#--coding:utf-8--
"""
callCuffquant.py
"""

__author__ = "CAO Yaqiang"
__date__ = "2019-07-01"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#sys
import os, time
from glob import glob
from datetime import datetime

#3rd
import click
from joblib import Parallel, delayed

#seqFlow
from utils import isTool, getLogger, callSys

#global setting
#logger
date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
logger = getLogger(fn=os.getcwd() + "/" + date.strip() + "_" +
                   os.path.basename(__file__) + ".log")


def callCuffquant(bam,gtf,cpu=3):
    fo = bam.split("/")[-2]
    if os.path.exists(fo):
        logger.info("%s has been generated to %s"%(bam,fo))
    cmd = "cuffquant -p {cpu} -o {fout} {gtf} {bam}".format(cpu=cpu,fout=fo,gtf=gtf,bam=bam)
    callSys([cmd],logger)


@click.command()
@click.option(
    "-pattern",
    required=True,
    help="Directory and patterns for the .bam files, for example './mouse*.bam'"
)
@click.option("-org",
              required=True,
              help="Organism for the data.",
              type=click.Choice(["hg38", "mm10"]))
@click.option("-cpu",
              default=10,
              help="Number of CPUs to finish the job, default is set to 10.")
def main(pattern,cpu=5,org="mm10"):
    for t in ["cuffquant"]:
        if not isTool(t):
            logger.error("%s not exits!" % t)
            return
    if org == "hg38":
        GTF="/home/caoy7/caoy7/Projects/0.Reference/1.hg38/2.annotations/gencode_human_v30_exons.gtf"
    elif org == "mm10":
        GTF="/home/caoy7/caoy7/Projects/0.Reference/2.mm10/2.annotations/gencode_v21_exons.gtf"
    else:
        return
    fs = glob(pattern)
    cpu = min(cpu, len(fs))
    Parallel(n_jobs=cpu)(delayed(callCuffquant)(f,GTF) for f in fs)


if __name__ == '__main__':
    start_time = datetime.now()
    main()
    elapsed = datetime.now() - start_time
    print("The process is done")
    print("Time used:", elapsed)
