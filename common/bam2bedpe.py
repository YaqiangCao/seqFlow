#!/usr/bin/env python2.7
#--coding:utf-8--
"""
bam2bedpe.py
2019-05-23: basically finished
2019-05-28: updated as bedpe stats
2019-06-12: updated as keep low mapq records.
"""

__author__ = "CAO Yaqiang"
__date__ = "2019-05-23"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import os, time, commands, gzip
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


def bam2Bedpe(bam, bedpe):
    fd = os.path.splitext(bedpe)[0]
    d = os.path.dirname(bedpe)
    if not os.path.exists(d):
        os.mkdir(d)
    nb = bam.split("/")[-1]
    tmpbam = fd + ".2.bam"
    #important for paired end reads!!
    samsort = "samtools sort -n -@ 2 {bam} -T {pre} -o {tmpbam}".format(
        bam=bam, tmpbam=nb, pre=nb.replace(".bam", ""))
    rmunmaped = "samtools view -b -q 0 -F 4 {} >> {}".format(nb, tmpbam)
    callSys([samsort, rmunmaped], logger)
    bam2bedpe = "bamToBed -bedpe -i {bam} > {bedpe}".format(bam=tmpbam,
                                                            bedpe=bedpe)
    logger.info(bam2bedpe)
    status, output = commands.getstatusoutput(bam2bedpe)
    rmbam = "rm {} {}".format(tmpbam, nb)
    callSys([rmbam, "gzip %s" % bedpe], logger)


def main():
    """
    Batch converting from bam to bedpe.
    """
    bams = glob("../3.mapping/*/*.bam")
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
    print "The process is done"
    print "Time used:", elapsed
