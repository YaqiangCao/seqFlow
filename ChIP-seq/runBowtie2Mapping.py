#!/usr/bin/env python2.7
#--coding:utf-8--
"""
bowtie2Mapping.py
2018-08-09
"""

__author__ = "CAO Yaqiang"
__date__ = "2014-09-23"
__modified__ = "2015-03-20"
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import glob, os, time, commands
from datetime import datetime

#3rd library
import matplotlib as mpl
mpl.use("pdf")
mpl.rcParams["pdf.fonttype"] = 42
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import pylab
import brewer2mpl
import seaborn as sns

#my own
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


def prepare_fastq(Fastq_Root="2.Fastq/", ):
    """
    If the fastq files are well prepared, suitable. 
    """
    fastqs = glob.glob(Fastq_Root + "*.fastq")
    data = {}
    for fq in fastqs:
        s = os.path.split(fq)[1]
        s = s.replace(".fastq", "")
        if s.endswith("_1"):
            sample = s.replace("_1", "")
            if sample not in data:
                data[sample] = [0, 0]
            data[sample][0] = fq
        if s.endswith("_2"):
            sample = s.replace("_2", "")
            if sample not in data:
                data[sample] = [0, 0]
            data[sample][1] = fq
        if not s.endswith("_1") and not s.endswith("_2"):
            data[s] = [fq]
    return data


def bowtie_mapping(sample, fqs):
    #hg19Ind = "/picb/molsysbio/usr/caoyaqiang/1.Projects/2.WL_hESC/0.GenomeReference/1.hg19/3.BowtieTophat/1.Bowtie1/hg19"
    hg38Ind = "/picb/molsysbio/usr/caoyaqiang/1.Projects/17.Alu/1.GenomeReference/1.hg38/3.BowtieTophatIndex/2.Bowtie2Tophat2/1.Bowtie2/hg38"
    logger.info("Start mapping %s.\n" % sample)
    if not os.path.exists(sample):
        os.mkdir(sample)
    sam = sample + "/" + sample + ".sam"
    bam = sample + "/" + sample + ".bam"
    if len(fqs) == 1:
        #doBowtie="/picb/molsysbio/usr/caoyaqiang/1.Projects/17.Alu/0.Tools/bowtie-1.1.0/bowtie -p 12 -q --phred33-quals -y -k 1 -m 1 --chunkmbs 200 -S --best --seed 123 -shmem %s %s %s"%( hg38Ind,fqs[ 0 ],sam )
        doBowtie = "/picb/molsysbio/usr/caoyaqiang/1.Projects/17.Alu/0.Tools/bowtie2-2.2.3/bowtie2 -q -N 1 --end-to-end -k 1 -X 1e9 --no-mixed --no-discordant --no-dovetail --no-contain --no-overlap -x %s %s %s" % (
            hg38Ind, fqs[0], sam)
    else:
        #doBowtie="/picb/molsysbio/usr/caoyaqiang/1.Projects/17.Alu/0.Tools/bowtie-1.1.0/bowtie -p 12 -q --phred33-quals -y -k 1 -m 1 --chunkmbs 200 -S --best --seed 123 -shmem %s -1 %s -2 %s %s"%( hg38Ind,fqs[ 0 ],fqs[ 1 ], sam )
        doBowtie = "/picb/molsysbio/usr/caoyaqiang/1.Projects/17.Alu/0.Tools/bowtie2-2.2.3/bowtie2 -q -N 1 --end-to-end -k 1 -X 1e9 --no-mixed --no-discordant --no-dovetail --no-contain --no-overlap -x %s -1 %s -2 %s %s" % (
            hg38Ind, fqs[0], fqs[1], sam)
    samview = "samtools view -S %s -b -o %s" % (sam, bam)
    samsort = "samtools sort -@ 2 %s %s/%s" % (bam, sample, sample)
    samindex = "samtools index %s %s/%s.bai  " % (bam, sample, sample)
    rmsam = "rm %s" % (sam)
    cmds = [samview, samsort, samindex, rmsam]
    logger.info("FLAG_A " + doBowtie)
    status, output = commands.getstatusoutput(doBowtie)
    #trim with "Warning"
    output = output.split("\n")
    output = [t for t in output if not t.startswith("Warning")]
    output = "\n".join(output)
    logger.info(sample + ":\n" + output + "\n")
    call_sys(cmds)
    return bam


def rm_dup_unmap(bam):
    fd = os.path.splitext(bam)[0]
    sam = fd + ".sam"
    bai = fd + ".bai"
    bed = fd + ".bed"
    tdf = fd + ".tdf"
    tmpbam = fd + "2.bam"
    tmpbai = fd + "2.bai"
    #get head
    cmd = "samtools view -H {}".format(bam)
    status, head = commands.getstatusoutput(cmd)
    with open(sam, "w") as f:
        f.write(head + "\n")
    #remove pcr duplicates
    rmdup = "samtools rmdup -s {} {}".format(bam, tmpbam)
    samindex1 = "samtools index {} {}".format(tmpbam, tmpbai)
    #remove unmaped reads
    rmunmaped = "samtools view -F 4 {} >> {}".format(tmpbam, sam)
    #convert sam to bam and index
    samview = "samtools view -S {} -b -o {}".format(sam, bam)
    samindex2 = "samtools index {} {}".format(bam, bai)
    rm = "rm {} {} {}".format(sam, tmpbam, tmpbai)
    #convert bam to bed, keep no-redundance
    bam2bed = "bamToBed -i {} > {}".format(bam, bed)
    #cmds = [ rmdup,samindex1,rmunmaped,samview,samindex2,rm,bam2bed ]
    cmds = [rmdup, samindex1, rmunmaped, samview, samindex2, rm]
    call_sys(cmds)
    status, output = commands.getstatusoutput("samtools flagstat %s" % bam)
    logger.info(
        "%s has been removed duplicates and unmapped reads, basic statistic: \n%s"
        % (bam, "FLAG_B " + "\n".join(output.split("\n")[:3])))


def trimBed(bed):
    #trim the 4th column of raw bed file to thin it
    n = bed + ".2"
    with open(n, "w") as f:
        for line in open(bed):
            line = line.split("\n")[0].split("\t")
            try:
                line[3] = line[5]
            except:
                report = "FLAG_C an error line in %s: %s" (bed,
                                                           "\t".join(line))
                continue
            line = "\t".join(line) + "\n"
            f.write(line)
    os.system("mv %s %s" % (n, bed))


def mapping(sample, fqs):
    bed = sample + ".bed"
    bam = sample + ".bam"
    if os.path.exists(sample + "/" + bam):
        #print bam,"already generated."
        return
    bam = bowtie_mapping(sample, fqs)
    rm_dup_unmap(bam)
    """
    #decompress the fq.gz
    cmd = "gunzip -c %s > %s"%( fq, nfq)
    call_sys( [ cmd ] )
    cmd = "rm %s"%nfq
    call_sys( [ cmd ] )
    trimBed( bed )
    fd = bam.split( "/" )[ 0 ]
    cmd1 = "mv %s ./"%bed
    cmd2 = "rm -r %s"%fd
    cmd3 = "bzip2 -z %s"%( fn+".bed" )
    cmd4 = "rm %s"%nfq
    cmds = [ cmd1, cmd2,cmd3,cmd4 ]
    for c in cmds:
        print c
        os.system( c )
    """


def bgzip_tabix(bedbz2):
    """
    Convert previously compressed bz2 file (bed) into bzip compressed file and index by tabix.
    """
    bed = bedbz2.replace(".bz2", "")
    bedgz = bed + ".gz"
    tbi = bedgz + ".tbi"
    if os.path.exists(bedgz) and os.path.exists(tbi):
        print bedgz, tbi, "has beed generated."
        return
    c1 = "bzip2 -d %s" % bedbz2
    c2 = "bgzip %s" % bed
    c3 = "tabix -s 1 -b 2 -e 3 %s" % bedgz
    call_sys([c1, c2, c3])


def sParseBowtie(lines):
    d, s = None, None
    lines = lines.split("\n")[-6:]
    s = lines[0].split()[-1][:-1]
    totalReads = int(lines[1].split(":")[1])
    d1 = lines[2].split(":")[1]
    uniqueMappedReads = int(d1.split("(")[0])
    uniqueMappedRatio = float(d1.split("(")[1][:-2])
    d2 = lines[4].split(":")[1]
    multipleMappedReads = int(d2.split("(")[0])
    multipleMappedRatio = float(d2.split("(")[1][:-2])
    d = {
        "TotalReads": totalReads,
        "UniqueMapReads": uniqueMappedReads,
        "UniqueMapRatio": uniqueMappedRatio,
        "MultipleMapReads": multipleMappedReads,
        "MultipleMapRatio": multipleMappedRatio,
    }
    #print s, d
    return d, s


def parseBowtielog(logs=None):
    if logs == None:
        logs = glob.glob("*.log")
    data = {}
    for log in logs:
        lines = open(log).read().split("\n\n")
        lines = [line for line in lines if "#" in line]
        for line in lines:
            if "#" in line:
                d, s = sParseBowtie(line)
                data[s] = d
    data = pd.DataFrame(data).T
    return data


def plotStat(data, pre="MappingStat"):
    fig, ax = pylab.subplots(figsize=(30, 6))
    colors = brewer2mpl.get_map('Set3', 'qualitative', 8).mpl_colors
    width = 0.5
    samples = list(data.index)
    ind = np.arange(len(data.index))
    #total reads as bar plot, seperate polyA and NonPolyA
    rs = data.loc[:, "TotalReads"].values
    rs = rs / 1000000.0
    ax.bar(ind, rs, width, color=colors[0])
    ax.set_ylabel("Reads/Million")
    ax.set_xticks(ind + width)
    ax.set_xticklabels(samples, rotation=45, ha="right")
    #map ratio
    ax2 = ax.twinx()
    umr = data.loc[:, "UniqueMapRatio"].values
    mmr = data.loc[:, "MultipleMapRatio"].values
    ax2.plot(
        ind + width / 2,
        umr,
        color="b",
        marker="s",
        label="Unique Mapping Ratio ")
    ax2.plot(
        ind + width / 2,
        mmr,
        color="r",
        marker="h",
        label="Multiple Mapping Ratio")
    for t in ax2.get_yticklabels():
        t.set_color("r")
    ax2.set_ylim([0, 100])
    ax2.set_ylabel("Maping Ratio/%")
    leg = ax2.legend(loc="upper right", fancybox=True)
    leg.get_frame().set_alpha(0.5)
    pylab.savefig(pre + ".pdf", dpi=1000, bbox_inches="tight")


def main():
    data = prepare_fastq(Fastq_Root="../1.RemoveAdapters/")
    Parallel(n_jobs=5)(delayed(mapping)(sample, fqs)
                       for sample, fqs in data.items())
    data = parseBowtielog()
    data.to_csv("MappingStat.txt", sep="\t", index_label="samples")
    plotStat(data)


if __name__ == '__main__':
    start_time = datetime.now()
    main()
    elapsed = datetime.now() - start_time
    print "The process is done"
    print "Time used:", elapsed
