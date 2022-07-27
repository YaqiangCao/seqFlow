#!/usr/bin/env python
#--coding:utf-8--
"""
CIGAR string https://www.drive5.com/usearch/manual/cigar.html
Very important for RNA reads
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


def get(name, key="RNA"):
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


def doMapping(sample, fq, starIndex, cpus=25):
    doStar = "STAR --runMode alignReads --runThreadN {cpus} --genomeDir {genomeDir} --readFilesIn {fq} --readFilesCommand 'zcat -1' --outFileNamePrefix {prefix}_ --outSAMtype BAM Unsorted --outSAMstrandField intronMotif --outFilterMultimapNmax 1".format(
        genomeDir=starIndex, cpus=cpus, fq=fq, prefix=sample)
    callSys([doStar], logger)
    #test for local
    c = "rm -fvr {pre}_STARtmp".format(pre=sample)
    #sorting and index bam, bzip unmapped fastq files
    bam = sample + "_Aligned.out.bam"
    bai = sample + "_Aligned.out.bai"
    pre = sample + "_Aligned.out"
    samsort = "samtools sort {bam} -T {pre} -o {bam}".format(bam=bam, pre=pre)
    samindex = "samtools index {bam} {bai}".format(bam=bam, bai=bai)
    cmds = [samsort, samindex]
    callSys(cmds, logger)


def bam2bed(sample, bam):
    """
    SAM to BAM and bedpe file.
    """
    bam2bed = "bamToBed -i {bam} -cigar > {bed}".format(bam=bam,
                                                        bed=sample +
                                                        "_RNA.bed")
    cmd = "gzip {bed}".format(bed=sample + "_RNA.bed")
    callSys([bam2bed, cmd], logger)


@click.command()
@click.option("-name", required=True, help="Sample name/id for the data.")
@click.option("-org",
              required=True,
              help="Organism for the data.",
              type=click.Choice(["hg38", "mm10"]))
def main(name, org):
    get(name, key="RNA")
    fq = "%s_RNA_R1.fastq.gz" % name
    starRef = {
        "hg38":
        "/home/caoy7/caoy7/Projects/0.Reference/1.hg38/3.index/5.star/2.star_2.7",
        "mm10":
        "/home/caoy7/caoy7/Projects/0.Reference/2.mm10/3.index/3.star/2.star_2.7",
    }
    sizeRef = {
        "hg38":
        "/home/caoy7/caoy7/Projects/0.Reference/1.hg38/1.fa/hg38.chrom.sizes",
        "mm10":
        "/home/caoy7/caoy7/Projects/0.Reference/2.mm10/1.fa/mm10.chrom.sizes",
    }
    ref = starRef[org]
    size = sizeRef[org]
    doMapping(name, fq, ref, cpus=25)
    bam2bed(name, name + "_Aligned.out.bam")


if __name__ == '__main__':
    main()
