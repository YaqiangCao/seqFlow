#!/usr/bin/env python2.7
#--coding:utf-8--
"""
tracMapping.py
2019-05-20
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

#my own
from utils import getLogger,callSys

#global settings
#logger
date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
logger = getLogger(fn=os.getcwd() + "/" + date.strip() + "_" +
                   os.path.basename(__file__) + ".log")




def prepare_fastq(Fastq_Root="../2.reid/"):
    """
    If the fastq files are well prepared, suitable. 
    """
    fastqs = glob(Fastq_Root + "*.fastq.gz")
    data = {}
    for fq in fastqs:
        s = os.path.split(fq)[1]
        s = s.replace(".fastq.gz", "")
        if s.endswith("_R1"):
            sample = s.replace("_R1", "")
            if sample not in data:
                data[sample] = [0, 0]
            data[sample][0] = fq
        if s.endswith("_R2"):
            sample = s.replace("_R2", "")
            if sample not in data:
                data[sample] = [0, 0]
            data[sample][1] = fq
        if not s.endswith("_R1") and not s.endswith("_R2"):
            data[s] = [fq]
    return data


def sam2bam(sam, bam):
    """
    SAM to BAM file 
    """
    samview = "samtools view -S %s -b -o %s" % (sam, bam)
    samsort = "samtools sort -@ 2 {bam} -T {pre} -o {bam}".format(
        bam=bam, pre=bam.replace(".bam", ""))
    samindex = "samtools index {bam} {bai}".format(bam=bam,
                                                   bai=bam.replace(
                                                       ".bam", ".bai"))
    rmsam = "rm %s" % (sam)
    cmds = [samview, samsort, samindex, rmsam]
    callSys(cmds,logger)


def tracMapping(sample, fqs, ref,cpus=10):
    logger.info("Start mapping %s.\n" % sample)
    if not os.path.exists(sample):
        os.mkdir(sample)
    sam = sample + "/" + sample + ".sam"
    bam = sample + "/" + sample + ".bam"
    if os.path.isfile(bam):
        logger.info("%s:%s exists! return." % (sample, bam))
        return
    if len(fqs) == 1:
        doBowtie = "bowtie2 -p {cpus} -q -N 1 --local --very-sensitive -k 1 --no-devetail --no-contain --no-overlap -x {ref} {fq} -S {sam}".format(cpus=cpus,ref=ref,fq=fqs[0],sam=sam)
    else:
        doBowtie = "bowtie2 -p {cpus} -q -N 1 --local --very-sensitive -k 1 --no-dovetail --no-contain --no-overlap -x {ref} -1 {fq1} -2 {fq2} -S {sam}".format(cpus=cpus,ref=ref,fq1=fqs[0],fq2=fqs[1],sam=sam)
    logger.info(doBowtie)
    status, output = commands.getstatusoutput(doBowtie)
    #trim with "Warning"
    output = output.split("\n")
    output = [t for t in output if not t.startswith("Warning")]
    output = "\n".join(output)
    logger.info("FLAG_A:" + sample + "\n" + output + "\nFLAG_A\n")
    sam2bam(sam, bam)
    return bam


def sParseBowtie(lines):
    """
    Parse Bowtie2 log file, to obtain mapping stastics.
    """
    d, s = None, None
    lines = lines.split("\n")
    s = lines[0]
    totalReads = int(lines[1].split(";")[0].split()[0])
    d1 = lines[4].strip().split()
    conUniqueMappedReads = int(d1[0])
    d2 = lines[8].strip().split()
    unconUniqueMappedReads = int(d2[0])
    mapRatio = float(lines[15].split("%")[0])
    d = {
        "TotalReads": totalReads,
        "ConcordantlyUniqueMapReads": conUniqueMappedReads,
        "DisconcordantlyUniqueMapReads": unconUniqueMappedReads,
        "MappingRatio(%s)": mapRatio
        #"MultipleMapReads": multipleMappedReads,
        #"MultipleMapRatio": multipleMappedRatio,
    }
    return d, s


def parseBowtielog(logs=None):
    if logs == None:
        logs = glob("*.log")
    data = {}
    for log in logs:
        lines = open(log).read().split("FLAG_A\n")
        lines = [line for line in lines if "FLAG_A" in line]
        for line in lines:
            t = line.split("FLAG_A:")[1]
            d, s = sParseBowtie(t)
            data[s] = d
    data = pd.DataFrame(data).T
    return data



def main():
    data = prepare_fastq("../../2.reid/")
    ref = "/home/tangq/tangq/Projects/0.Reference/1.hg38/3.index/2.bowtie2/hg38"
    #Parallel(n_jobs=5)(delayed(tracMapping)(sample, fqs, ref)
    #                   for sample, fqs in data.items())
    data = parseBowtielog()
    data.to_csv("MappingStat.txt", sep="\t", index_label="samples")


if __name__ == '__main__':
    start_time = datetime.now()
    date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
    main()
    elapsed = datetime.now() - start_time
    print("The process is done")
    print("Time used:", elapsed)
