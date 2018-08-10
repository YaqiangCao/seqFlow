#!/usr/bin/env python2.7
#--coding:utf-8--
"""
callEditing.py
"""

__author__ = "CAO Yaqiang"
__date__ = "2018-08-08"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import os, gzip, time
from glob import glob

#3rd library
from joblib import Parallel, delayed

#my own
from Biolibs.rel.General.logger import getlogger

#global settings
#logger
date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
logger = getlogger(fn=os.getcwd() + "/" + date.strip() + "_" +
                   os.path.basename(__file__) + ".log")

#global 
giremi = "/picb/molsysbio2/caoyaqiang/2.WWX_RNA-seq/0.Tools/giremi-0.2.1/giremi"
FA = "/home/caoyaqiang/cyq_m2/2.WWX_RNA-seq/1.Reference/1.FAs/hg38.fa"


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


def preds():
    bamdir = "../../2.Mapping/"
    giremidir = "../1.giremi_l/"
    bams = glob(bamdir + "*/*.bam")
    gs = glob(giremidir + "*.giremi")
    ds = {}
    for bam in bams:
        n = bam.split("/")[-2]
        ds[n] = {"bam": bam}
    for g in gs:
        n = g.split("/")[-1].split(".giremi")[0]
        ds[n]["giremi"] = g
    for key in ds.keys():
        if "giremi" not in ds[key]:
            del ds[key]
    return ds


def runGerimi(key, fl, bam):
    fout = key + ".giremi"
    tmp = key + ".tmp"
    if os.path.isfile(tmp) or os.path.isfile(fout):
        return
    with open(tmp, "w") as fo:
        pass
    cmd1 = "{giremi} -f {fa} -o {fout} -l {giremi_l} -m 20 -s 0 {bam}".format(
        giremi=giremi, fa=FA, fout=fout, giremi_l=fl, bam=bam)
    cmd2 = "rm %s" % tmp
    call_sys([cmd1, cmd2])


def main():
    ds = preds()
    Parallel(n_jobs=10)(delayed(runGerimi)(key,
                                           ds[key]["giremi"], ds[key]["bam"])
                        for key in ds.keys())


main()
