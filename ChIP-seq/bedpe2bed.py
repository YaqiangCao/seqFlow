#!/usr/bin/env python2.7
#--coding:utf-8 --
"""
bedpe2bed.py
"""

__author__ = "CAO Yaqiang"
__date__ = "2015-05-26"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#general library
import glob, os, time, sys, random

#3rd library

#my own library
from Biolibs.rel.General.logger import getlogger

#global settings
#logger
date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
logger = getlogger(fn=os.getcwd() + "/" + date.strip() + "_" +
                   os.path.basename(__file__) + ".log")


def call_sys(cmds):
    """
    Call systematic commands without return.
    """
    for c in cmds:
        logger.info(c)
        try:
            os.system(c)
        except:
            logger.error(c)


def flush(i):
    if i % 1000 == 0:
        report = "\r%dk reads parsed" % (i / 1000)
        sys.stdout.write(report)
        sys.stdout.flush()


def bedpe2bed(bedpe):
    fn = os.path.splitext(bedpe)[0]
    bed = fn + '.bed'
    with open(bed, "w") as f:
        for i, line in enumerate(open(bedpe)):
            flush(i)
            line = line.split("\n")[0].split("\t")
            if len(line) < 7:
                continue
            r1 = [line[0], line[1], line[2], line[6]]
            r2 = [line[3], line[4], line[5], line[6]]
            r1 = "\t".join(r1) + "\n"
            r2 = "\t".join(r2) + "\n"
            f.write(r1 + r2)
    print
    tmp = str(random.random())
    bedgz = bed + ".gz"
    c1 = "sortBed -i %s > %s" % (bed, tmp)
    c2 = "mv %s %s" % (tmp, bed)
    c3 = "bgzip %s" % bed
    c4 = "tabix -p bed %s" % bedgz
    call_sys([c1, c2, c3, c4])


def main():
    bedpe = "H3K27ac_chr22.bedpe"
    bedpe2bed(bedpe)


if __name__ == "__main__":
    main()
