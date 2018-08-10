#!/usr/bin/env python
#--coding:utf-8 --
"""
bamPileup.py
"""

import os, sys, shutil, gzip, glob, time, commands
from datetime import datetime
from joblib import Parallel, delayed

from Biolibs.rel.General.logger import getlogger
date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
logger = getlogger(fn=os.getcwd() + "/" + date.strip() + "_" +
                   os.path.basename(__file__) + ".log")

CHROM_SIZE = "/picb/molsysbio/usr/caoyaqiang/1.Projects/5.hESCSortedlncRNAsPipeline/1.GenomeReference/1.Fas/hg19.chrom.sizes"

__author__ = "CAO Yaqiang"
__date__ = "2014-09-17"
__modified__ = "2015-01-05"
__email__ = "caoyaqiang0410@gmail.com"


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


def get_chr_size():
    chrs = {}
    for line in open(CHROM_SIZE):
        line = line.split("\n")[0].split("\t")
        chrs[line[0]] = int(line[1])
    return chrs


def validate_bdg(bdg):
    chrs = get_chr_size()
    nbdg = bdg + ".2"
    with open(nbdg, "w") as f:
        for line in open(bdg):
            line = line.split("\n")[0].split("\t")
            if len(line) < 4:
                continue
            if int(line[1]) >= chrs[line[0]]:
                line[1] = str(chrs[line[0]] - 1)
            if int(line[2]) >= chrs[line[0]]:
                line[2] = str(chrs[line[0]] - 1)
            line = "\t".join(line) + "\n"
            f.write(line)
    cmd = "mv %s %s" % (nbdg, bdg)
    os.system(cmd)


def bdg2bw(bdg):
    print "Validating :", bdg
    validate_bdg(bdg)
    bw = os.path.splitext(bdg)[0] + ".bw"
    cmd1 = "bedGraphToBigWig {bamgraph} {chromsize} {bigwig}".format(
        bamgraph=bdg, chromsize=CHROM_SIZE, bigwig=bw)
    cmd2 = "rm %s" % bdg
    call_sys([cmd1, cmd2])


def after_macs2_pileup(dir):
    f = glob.glob(dir + "/*.r")[0]
    cmdr = "Rscript %s" % f
    f = glob.glob(dir + "/*control*.bdg")[0]
    cmdrm = "rm %s" % f
    call_sys([cmdr, cmdrm])
    bdg = glob.glob(dir + "/*treat*.bdg")[0]
    nbdg = bdg.split("_treat")[0] + ".bdg"
    cmd = "mv %s %s" % (bdg, nbdg)
    call_sys([cmd])


def raw_pileup(bam, sample, root_dir="1.RAW_Pileup/"):
    dir = os.path.join(root_dir, sample)
    cmdmacs2 = "macs2 callpeak -t {} -n {} -B --outdir {} --fix-bimodal".format(
        bam, sample, dir)
    call_sys([cmdmacs2])
    after_macs2_pileup(dir)


def rpm_pileup(bam, sample, root_dir="2.RPM_Pileup/"):
    dir = os.path.join(root_dir, sample)
    cmdmacs2 = "macs2 callpeak -t {} -n {} -B --SPMR --outdir {} --fix-bimodal".format(
        bam, sample, dir)
    call_sys([cmdmacs2])
    after_macs2_pileup(dir)


def pileup(bam, root_raw, root_rpm):
    sample = os.path.splitext(os.path.split(bam)[1])[0]
    raw_bdg = os.path.join(root_raw, sample, sample + ".bdg")
    rpm_bw = os.path.join(root_rpm, sample, sample + ".bw")
    if os.path.exists(raw_bdg) and os.path.exists(rpm_bw):
        print raw_bdg, rpm_bw, "has been generated"
        return
    if not os.path.exists(raw_bdg):
        try:
            raw_pileup(bam, sample, root_dir=root_raw)
        except:
            logger.info("ERROR!%s" % sample)
    if not os.path.exists(rpm_bw):
        try:
            rpm_pileup(bam, sample, root_dir=root_rpm)
        except:
            logger.info("ERROR!%s" % sample)


def prepareBdgs(root_dir="1.RAW_Pileup"):
    bdgs = glob.glob(os.path.join(root_dir, "*/*.bdg"))
    samples = {}
    for bdg in bdgs:
        n = os.path.splitext(os.path.split(bdg)[1])[0]
        treat = n.split("_")[0]
        antibody = n.split("_")[1]
        if treat not in samples:
            samples[treat] = {"input": None, "antibody": []}
        if antibody == "input":
            samples[treat]["input"] = bdg
        else:
            samples[treat]["antibody"].append(bdg)
    nbdgs = []
    for s in samples:
        for t in samples[s]["antibody"]:
            nbdgs.append((s, t, samples[s]["input"]))
    return nbdgs


def trackFE(sample, treatbdg, inputbdg, root_dir="3.FE_Pileup/"):
    fn = os.path.split(treatbdg)[-1]
    outbdg = os.path.join(root_dir, fn)
    cmd = "macs2 bdgcmp -t {treatbdg} -c {inputbdg} -o {outbdg} -m logFE -p 1e-6".format(
        treatbdg=treatbdg, inputbdg=inputbdg, outbdg=outbdg)
    call_sys([cmd])


def main():
    mapping_root = "../3.Mapping"
    root_raw, root_rpm, root_fe = "1.RAW_Pileup", "2.RPM_Pileup", "3.FE_Pileup"
    rs = [root_raw, root_rpm, root_fe]
    for r in rs:
        if not os.path.exists(r):
            try:
                os.mkdir(r)
            except:
                pass
    #bams = glob.glob( os.path.join(mapping_root,"*/*.bam") )
    #Parallel( n_jobs=30 )( delayed( pileup )( f, root_raw, root_rpm) for f in bams )
    #bdgs = prepareBdgs(  )
    #Parallel( n_jobs=28 )( delayed( trackFE )( s,t,i,root_fe ) for s,t,i in bdgs )
    bdgs = glob.glob("1.RAW_Pileup/*/*.bdg")
    bdgs.extend(glob.glob("2.RPM_Pileup/*/*.bdg"))
    bdgs.extend(glob.glob("3.*/*.bdg"))
    Parallel(n_jobs=30)(delayed(bdg2bw)(b) for b in bdgs)


if __name__ == "__main__":
    startTime = datetime.now()
    main()
    elapsed = datetime.now() - startTime
    print "The process is done."
    print "Time used:", elapsed
