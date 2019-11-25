#!/usr/bin/env python
#--coding:utf-8--
"""
bdg2bw.py
2019-06-24: updated
"""

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

CHROM = None


def getChrSize():
    chrs = {}
    for line in open(CHROM):
        line = line.split("\n")[0].split("\t")
        chrs[line[0]] = int(line[1])
    return chrs


def validateBdg(bdg):
    """
    Validating .bdg files through chrom size.
    """
    chrs = getChrSize()
    nbdg = bdg + ".2"
    with open(nbdg, "w") as f:
        for line in open(bdg):
            line = line.split("\n")[0].split("\t")
            if len(line) < 4:
                continue
            if line[0] not in chrs:
                continue
            if int(line[1]) >= chrs[line[0]] or int(line[2]) > chrs[line[0]]:
                continue
            line = "\t".join(line) + "\n"
            f.write(line)
    cmd = "mv %s %s" % (nbdg, bdg)
    callSys([cmd], logger)


def bdg2bw(f):
    """
    Converting .bdg file to .bw file through bedGraphToBigWig.
    """
    n = f.split("/")[-1].replace(".bdg", "")
    if os.path.isfile(n + ".bw"):
        return
    cmd1 = "bedSort {bdg} {sbdg}".format(bdg=f, sbdg=n + ".bdg2")
    callSys([cmd1], logger)
    #validation bdg
    validateBdg(n + ".bdg2")
    cmd2 = "bedGraphToBigWig {bdg} {chrom} {bw}".format(bdg=n + ".bdg2",
                                                        chrom=CHROM,
                                                        bw=n + ".bw")
    callSys([cmd2, "rm %s.bdg2" % n], logger)


@click.command()
@click.option(
    "-pattern",
    required=True,
    help="Directory and patterns for the .bg files, for example './mouse*.bdg'"
)
@click.option("-org",
              required=True,
              help="Organism for the data.",
              type=click.Choice(["hg38", "mm10"]))
@click.option("-cpu",
              default=10,
              help="Number of CPUs to finish the job, default is set to 10.")
def main(pattern, org, cpu):
    global CHROM
    for t in ["bedSort", "bedGraphToBigWig"]:
        if not isTool(t):
            logger.error("%s not exits!" % t)
            return
    if org == "hg38":
        CHROM = "/home/caoy7/code/seqFlow/data/hg38.chrom.sizes"
    elif org == "mm10":
        CHROM = "/home/caoy7/code/seqFlow/data/mm10.chrom.sizes"
    else:
        return
    fs = glob(pattern)
    cpu = min(cpu, len(fs))
    Parallel(n_jobs=cpu)(delayed(bdg2bw)(f) for f in fs)


if __name__ == "__main__":
    startTime = datetime.now()
    main()
    elapsed = datetime.now() - startTime
