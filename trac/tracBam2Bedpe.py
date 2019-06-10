#!/usr/bin/env python2.7
#--coding:utf-8--
"""
tracPreBam.py
2019-05-23
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


def bam2Bedpe(bam, bedpe, mapq=10):
    fd = os.path.splitext(bedpe)[0]
    d = os.path.dirname(bedpe)
    if not os.path.exists(d):
        os.mkdir(d)
    tmpbam = fd + "2.bam"
    rmunmaped = "samtools view -q 10 -b -F 4 {} >> {}".format(bam, tmpbam)
    callSys([rmunmaped], logger)
    bam2bedpe = "bamToBed -bedpe -i {bam} > {bedpe}".format(bam=tmpbam,
                                                            bedpe=bedpe)
    logger.info(bam2bedpe)
    status, output = commands.getstatusoutput(bam2bedpe)
    rmbam = "rm {}".format(tmpbam)
    callSys([rmbam], logger)


def main():
    bams = glob("../3.mapping/*/*.bam")
    ds = []
    for bam in bams:
        nb = "./" + bam.split("/")[-1].replace(".bam", ".bedpe")
        if os.path.isfile(nb):
            logger.info("%s has been generated. return." % nb)
            continue
        ds.append([bam, nb])
    Parallel(n_jobs=5)(delayed(bam2Bedpe)(t[0], t[1]) for t in ds[:1])
    fs = glob("*.bedpe")
    Parallel(n_jobs=5)(delayed(thinBedpe)(f) for f in fs)


if __name__ == '__main__':
    start_time = datetime.now()
    main()
    elapsed = datetime.now() - start_time
    print "The process is done"
    print "Time used:", elapsed
