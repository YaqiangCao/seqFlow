#!/usr/bin/env python2.7
#--coding:utf-8 --
"""

"""

__author__ = "CAO Yaqiang"
__date__ = ""
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#general library
import os, time, sys, random
from glob import glob

#3rd library
from joblib import Parallel, delayed

#trac
from utils import getLogger, callSys, cFlush

#global settings
#logger
date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
logger = getLogger(fn=os.getcwd() + "/" + date.strip() + "_" +
                   os.path.basename(__file__) + ".log")


def bedpe2bed(bedpe):
    fn = os.path.splitext(os.path.split(bedpe)[1])[0]
    bed = fn + '.bed'
    with open(bed, "w") as f:
        for i, line in enumerate(open(bedpe)):
            cFlush(i)
            line = line.split("\n")[0].split("\t")
            if len(line) < 7:
                continue
            if "_" in line[0]:
                continue
            r1 = [line[0], line[1], line[2], line[6]]
            r2 = [line[3], line[4], line[5], line[6]]
            r1 = "\t".join(r1) + "\n"
            r2 = "\t".join(r2) + "\n"
            f.write(r1 + r2)
    tmp = str(random.random())
    bedgz = bed + ".gz"
    c1 = "sortBed -i %s > %s" % (bed, tmp)
    c2 = "mv %s %s" % (tmp, bed)
    c3 = "bgzip %s" % bed
    c4 = "tabix -p bed %s" % bedgz
    callSys([c1, c2, c3, c4], logger)


def main():
    fs = glob("../../4.uniqueNonRedudantBEDPE/2.bedpe/*.bedpe")
    data = Parallel(n_jobs=len(fs))(delayed(bedpe2bed)(f) for f in fs)


if __name__ == "__main__":
    main()
