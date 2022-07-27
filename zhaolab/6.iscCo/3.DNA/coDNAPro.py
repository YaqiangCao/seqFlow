#!/usr/bin/env python
#--coding:utf-8--
"""
"""

__author__ = "CAO Yaqiang"
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import os
import time
import subprocess
from glob import glob

#3rd library
import click
from joblib import Parallel, delayed

#pipi
from pipi.utils import getLogger, callSys, isTool

#global settings
#logger
date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
logger = getLogger(fn=os.getcwd() + "/" + date.strip() + "_" +
                   os.path.basename(__file__) + ".log")


def get(name, key="DNA"):
    fs = glob("../2.sepSingleDNARNA/*%s_R1.fastq.gz" % key)
    fs.sort()
    cmd = "cat %s > %s_%s_R1.fastq.gz" % (" ".join(fs), name, key)
    callSys([cmd], logger)


def doMapping(sample, fq, ref, outdir, cpus=25):
    """
    Mapping settings for iscDNase-seq data.
    """
    od = outdir
    sam = od + "/" + sample + ".sam"
    bam = od + "/" + sample + ".bam"
    #cpus=cpus, ref=ref, fq1=fqs[0], fq2=fqs[1], sam=sam)
    doBowtie = "bowtie2 --no-mixed --no-discordant -p {cpus} -q --local --very-sensitive -x {ref} {fq} -S {sam}".format(
        cpus=cpus, ref=ref, fq=fq, sam=sam)
    logger.info(doBowtie)
    stat, output = subprocess.getstatusoutput(doBowtie)
    #trim with "Warning"
    output = output.split("\n")
    output = [t for t in output if not t.startswith("Warning")]
    output = "\n".join(output)
    logger.info("FLAG_A:" + sample + "\n" + output + "\nFLAG_A\n")
    return sample, sam


def sam2bamBed(sample, sam, mapq=10):
    """
    SAM to BAM and bedpe file.
    """
    n = os.path.splitext(sam)[0]
    bam = n + ".bam"
    bedAll = n + ".bed"
    #sam to bam, filtering mapq
    samview = "samtools view -b -F 4 -@ 2 -q {mapq} -o {bam} {sam}".format(
        mapq=mapq, bam=bam, sam=sam)
    #sort by read name
    #samsort = "samtools sort -n -@ 2 {bam} -T {pre} -o {bam}".format(
    samsort = "samtools sort -@ 2 {bam} -T {pre} -o {bam}".format(
        bam=bam, pre=bam.replace(".bam", ""))
    rmsam = "rm %s" % (sam)
    cmds = [samview, samsort, rmsam]
    callSys(cmds, logger)
    bam2bed = "bamToBed -i {bam} > {bed}".format(bam=bam, bed=bedAll)
    logger.info(bam2bed)
    stat, output = subprocess.getstatusoutput(bam2bed)
    cmd = "gzip %s" % (bedAll)
    callSys([cmd], logger)
    cmd = "samtools index {bam} {bai}".format(bam=bam,
                                              bai=bam.replace(".bam", ".bai"))
    callSys([cmd], logger)
    cmd = "bamCoverage -b {bam} -o {bw} -p 10 --ignoreDuplicates --minMappingQuality 10 --normalizeUsing CPM".format(
        bam=bam, bw=bam.replace(".bam", ".bw"))
    callSys([cmd], logger)
    return bedAll + ".gz"


@click.command()
@click.option("-name", required=True, help="Sample name/id for the data.")
@click.option("-org",
              required=True,
              help="Organism for the data.",
              type=click.Choice(["hg38", "mm10"]))
def main(name, org):
    get(name, key="DNA")
    fq = "%s_DNA_R1.fastq.gz" % name
    refData = {
        "hg38":
        "/home/caoy7/caoy7/Projects/0.Reference/1.hg38/3.index/2.bowtie2/hg38",
        "mm10":
        "/home/caoy7/caoy7/Projects/0.Reference/2.mm10/3.index/2.bowtie2/mm10",
    }
    ref = refData[org]
    sample, sam = doMapping("%s_DNA" % name, fq, ref, "./", cpus=25)

    #step 3, convert to bam and bedpe files
    #sam to bam and bedpe
    bed = sam2bamBed("%s_DNA" % name, "%s_DNA.sam" % name)


if __name__ == '__main__':
    main()
