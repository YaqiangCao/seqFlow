#!/usr/bin/env python2.7
#--coding:utf-8--
"""
champBam2Bedpe.py
"""

__author__ = "CAO Yaqiang"
__date__ = "2019-12-13"
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import os, time, subprocess,gzip 
from glob import glob
from datetime import datetime
from collections import Counter

#3rd library
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

#this
#from utils import getLogger, callSys, PET
from utils import *

#global settings
#logger
date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
logger = getLogger(fn=os.getcwd() + "/" + date.strip() + "_" +
                   os.path.basename(__file__) + ".log")


def bam2Bedpe(bam, bedpe, mapq=10):
    """
    Converting BAM file to BED file. 
    bam: bam file path
    bedpe: bedpe file path
    mapq: mapq cutoff to remove bad qulity reads.
    """
    fd = os.path.splitext(bedpe)[0]
    d = os.path.dirname(bedpe)
    if not os.path.exists(d):
        os.mkdir(d)
    nb = bam.split("/")[-1]
    tmpbam = fd + ".2.bam"
    #important for paired end reads!!
    rmunmaped = "samtools view -b -q {} -F 4 {} >> {}".format(
        mapq, bam, tmpbam)
    samsort = "samtools sort -n -@ 2 {bam} -T {pre} -o {tmpbam}".format(
        bam=tmpbam, tmpbam=nb, pre=nb.replace(".bam", ""))
    callSys([rmunmaped, samsort], logger)
    bam2bedpe = "bamToBed -bedpe -i {bam} > {bedpe}".format(bam=nb,
                                                            bedpe=bedpe)
    logger.info(bam2bedpe)
    stat, output = subprocess.getstatusoutput(bam2bedpe)
    rmbam = "rm {} {}".format(tmpbam, nb)
    callSys([rmbam, "gzip %s" % bedpe], logger)


def main():
    """
    Batch converting from bam to bedpe.
    """
    bams = glob("../2.mapping/*/*.bam")
    ds = []
    for bam in bams:
        bai = bam.replace(".bam", ".bai")
        if not os.path.isfile(bai):
            continue
        nb = "./" + bam.split("/")[-1].replace(".bam", ".bedpe")
        if os.path.isfile(nb):
            logger.info("%s has been generated. return." % nb)
            continue
        if os.path.isfile(nb + ".gz"):
            logger.info("%s has been generated. return." % (nb + ".gz"))
            continue
        ds.append([bam, nb])
    Parallel(n_jobs=20)(delayed(bam2Bedpe)(t[0], t[1]) for t in ds)


if __name__ == '__main__':
    start_time = datetime.now()
    main()
    elapsed = datetime.now() - start_time
    print("The process is done")
    print("Time used:", elapsed)
