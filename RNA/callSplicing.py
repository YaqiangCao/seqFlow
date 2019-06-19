#!/usr/bin/env python2.7
#--coding:utf-8--
"""
callSplicing.py
2015-01-15: If packed, not fit for ploting.
2015-01-21: Modified DEIs filtering cutoffs, custome cutoff as psi range cutoff added. 
2015-01-23: Modified as added all alternative splicing types.
"""

__author__ = "CAO Yaqiang"
__date__ = "2015-01-23"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import glob, os, time, commands, shutil
from datetime import datetime

#3rd library
import pandas as pd
from joblib import Parallel, delayed

#my own
from Biolibs.rel.General.logger import getlogger

#global settings
#logger
date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
logger = getlogger(fn=os.getcwd() + "/" + date.strip() + "_" +
                   os.path.basename(__file__) + ".log")
#data
#index gff file
INDEXS = {
    "SE":
    "/picb/molsysbio/usr/caoyaqiang/1.Projects/17.Alu/1.GenomeReference/9.MISO/3.Index/SE",
    "RI":
    "/picb/molsysbio/usr/caoyaqiang/1.Projects/17.Alu/1.GenomeReference/9.MISO/3.Index/RI",
    "A3SS":
    "/picb/molsysbio/usr/caoyaqiang/1.Projects/17.Alu/1.GenomeReference/9.MISO/3.Index/A3SS",
    "A5SS":
    "/picb/molsysbio/usr/caoyaqiang/1.Projects/17.Alu/1.GenomeReference/9.MISO/3.Index/A5SS",
    "MXE":
    "/picb/molsysbio/usr/caoyaqiang/1.Projects/17.Alu/1.GenomeReference/9.MISO/3.Index/MXE",
}
#generate by exon_utils --get-const-exons ../2.hg38/SE.hg38.gff3 --min-exon-size 1000  --output-dir ./
EXONS = "/picb/molsysbio/usr/caoyaqiang/1.Projects/17.Alu/1.GenomeReference/9.MISO/4.hg38_Exons/SE.hg38.min_1000.const_exons.gff"


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


def parseIDL(f):
    d = open(f).read().split("\n")[0][1:]
    d = d.split(",")
    mean = float(d[0].split("=")[1])
    std = float(d[1].split("=")[1])
    return mean, std


def pipeline(bam, rlen=111, rs=[
        "1.ILD",
        "2.PSI",
]):
    #specific for STAR output bam files
    pre = os.path.split(bam)[1]
    pre = os.path.splitext(pre)[0]
    pre = pre.split("_Aligned")[0]
    ild_root = os.path.join(rs[0], pre)
    c1 = "pe_utils --compute-insert-len {bam} {exons} --output {out}".format(
        bam=bam, exons=EXONS, out=ild_root)
    call_sys([c1])
    f = glob.glob(ild_root + "/*.insert_len")[0]
    mean, std = parseIDL(f)
    #run each alternative splicing events
    for key, value in INDEXS.items():
        nRoot = os.path.join(rs[1], key, pre)
        c2 = "miso --run {index} {bam} --output-dir {out} --read-len {rlen} --paired-end {mean} {std} --prefilter -p 5".format(
            index=value, bam=bam, out=nRoot, rlen=rlen, mean=mean, std=std)
        c3 = "summarize_miso --summarize-samples {root} {root}".format(
            root=nRoot)
        #c4 = "miso_pack --pack {root}".format( root=se_root )
        call_sys([
            c2,
            c3,
        ])


def getDEIs(s1, s2, pre, rs=["3.DEIs", "4.Filter"]):
    pre1 = os.path.join(rs[0], pre)
    c1 = "compare_miso --compare-samples {s1} {s2} {pre}".format(s1=s1,
                                                                 s2=s2,
                                                                 pre=pre1)
    call_sys([c1])
    bf = glob.glob(pre1 + "/*/*/*.miso_bf")[0]
    pre2 = os.path.join(rs[1], pre)
    #the most strigent
    #c2 = "filter_events --filter {bf} --num-sum-inc-exc 20 --num-exc 3 --num-inc 3 --delta-psi 0.2 --bayes-factor 2 --apply-both --output-dir {pre}".format( bf=bf,pre=pre2 )
    #baysian factor set as 1, which is less strigent for following cutoff.
    c2 = "filter_events --filter {bf} --num-total 20 --delta-psi 0.2 --bayes-factor 1 --apply-both --output-dir {pre}".format(
        bf=bf, pre=pre2)
    call_sys([c2])


def furthurCut(deltaPsi=0.25, bf=1, root="4.Filter"):
    fs = glob.glob(root + "/*/*miso_bf.filtered")
    for f in fs:
        mat = pd.read_table(f, index_col=0)
        a, b, c, d = mat["sample1_ci_low"], mat["sample1_ci_high"], mat[
            "sample2_ci_low"], mat["sample2_ci_high"]
        g, h = b - a, d - c
        g = g[g < deltaPsi]
        h = h[h < deltaPsi]
        nis = set(g.index).intersection(set(h.index))
        e = mat["bayes_factor"]
        e = e[e > bf]
        nis = nis.intersection(set(e.index))
        mat = mat.loc[list(nis), :]
        mat.to_csv(f + "2", sep="\t", index_label="event")


def main():
    roots = ["1.ILD", "2.PSI", "3.DEIs", "4.Filter"]
    rlen = 86
    mapping_root = "../../3.Mapping"
    bams = glob.glob(mapping_root + "/*/*.bam")
    bams.sort()
    Parallel(n_jobs=6)(delayed(pipeline)(bam, rlen, roots[:2]) for bam in bams)


if __name__ == '__main__':
    start_time = datetime.now()
    main()
    elapsed = datetime.now() - start_time
    print "The process is done"
    print "Time used:", elapsed
