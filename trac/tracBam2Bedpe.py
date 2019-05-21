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


def thinBedpe(f):
    with open(f + ".2", "w") as fo:
        for i, line in enumerate(open(f)):
            line = line.split("\n")[0].split("\t")
            #shroten the name
            line[6] = str(i)
            fo.write("\t".join(line) + "\n")
    #cmds = ["mv %s %s" % (f + ".2", f), "gzip %s" % f]
    cmds = ["mv %s %s" % (f + ".2", f)]
    callSys(cmds, logger)


def bam2bedpe(bam):
    bedpe = bam.split("/")[-1].replace(".bam", ".bedpe")
    if os.path.isfile(bedpe):
        logger.info("%s has been generated! return" % bedpe)
        return
    cmd = "bamToBed -bedpe -i {bam} > {bedpe}".format(bam=bam, bedpe=bedpe)
    logger.info(cmd)
    status, output = commands.getstatusoutput(cmd)


def main():
    bams = glob("../1.bams/*.bam")
    Parallel(n_jobs=len(bams))(delayed(bam2bedpe)(f) for f in bams)
    fs = glob("*.bedpe")
    Parallel(n_jobs=len(fs))(delayed(thinBedpe)(f) for f in fs)


if __name__ == '__main__':
    start_time = datetime.now()
    main()
    elapsed = datetime.now() - start_time
    print "The process is done"
    print "Time used:", elapsed
