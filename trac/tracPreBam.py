#!/usr/bin/env python2.7
#--coding:utf-8--
"""
tracPreBam.py
"""

__author__ = "CAO Yaqiang"
__date__ = "2014-09-23"
__modified__ = "2015-03-20"
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import os, time, commands
from glob import glob
from datetime import datetime

#3rd library
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

#this
from utils import getLogger, callSys

#global settings
#logger
date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
logger = getLogger(fn=os.getcwd() + "/" + date.strip() + "_" +
                   os.path.basename(__file__) + ".log")


def getNodupMappedBam(bam, outbam):
    fd = os.path.splitext(outbam)[0]
    d = os.path.dirname(outbam)
    if not os.path.exists(d):
        os.mkdir(d)
    sam = fd + ".sam"
    bai = fd + ".bai"
    tmpbam = fd + "2.bam"
    tmpbai = fd + "2.bai"
    #get head
    cmd = "samtools view -H {}".format(bam)
    status, head = commands.getstatusoutput(cmd)
    with open(sam, "w") as f:
        f.write(head + "\n")
    #remove pcr duplicates
    rmdup = "samtools rmdup {} {}".format(bam, tmpbam)
    samindex1 = "samtools index {} {}".format(tmpbam, tmpbai)
    #remove unmaped reads
    rmunmaped = "samtools view -F 4 {} >> {}".format(tmpbam, sam)
    #convert sam to bam and index
    samview = "samtools view -S {} -b -o {}".format(sam, outbam)
    samindex2 = "samtools index {} {}".format(outbam, bai)
    rm = "rm {} {} {}".format(sam, tmpbam, tmpbai)
    cmds = [rmdup, samindex1, rmunmaped, samview, samindex2, rm]
    callSys(cmds, logger)
    status, output = commands.getstatusoutput("samtools flagstat %s" % bam)
    logger.info(
        "%s has been removed duplicates and unmapped reads, basic statistic: \n%s"
        % (bam, "FLAG_B " + "\n".join(output.split("\n")[:3])))


def main():
    bams = glob("../../3.mapping/1.local/*/*.bam")
    ds = []
    for bam in bams:
        nb = "./" + bam.split("/")[-1]
        if os.path.isfile(nb):
            logger.info("%s has been generated. return." % nb)
            continue
        ds.append([bam, nb])
    Parallel(n_jobs=5)(delayed(getNodupMappedBam)(t[0], t[1]) for t in ds)


if __name__ == '__main__':
    start_time = datetime.now()
    main()
    elapsed = datetime.now() - start_time
    print "The process is done"
    print "Time used:", elapsed
