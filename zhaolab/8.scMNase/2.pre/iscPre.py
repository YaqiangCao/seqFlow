#!/usr/bin/env python
#--coding:utf-8--
"""
iscPre.py 
Pre-processing indexing based split-seq single-cell reads to genomic reads.

Main functions:
1. Trim the barcode and linker, and rename read id 
2. Map to reference genome 
3. Process into unique TSV file for each read with cell barcode 
4. filter blacklist 
5. Show qc stats

2021-04-24: update parallel pre-processing for trim reads
2022-11-29: update unique reads counting, take cell id into consideration
2023-05-18: add blacklist filtering; using click as interface
"""

__author__ = "CAO Yaqiang"
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import os
import time
import gzip
import shutil
import random
import subprocess
from glob import glob
from datetime import datetime
from argparse import RawTextHelpFormatter

#3rd library
import click
import numpy as np
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator


def isTool(name):
    """
    Check if a tool is on PATH and marked as executable.

    @param name: str, 

    @return: bool, true or false
    """
    from distutils.spawn import find_executable
    if find_executable(name) is not None:
        return True
    else:
        return False


def runCmds(cs, p=True):
    """
    Run system commands
    """
    for c in cs:
        if p == True:
            print("\t".join([str(datetime.now()), c]))
        output = subprocess.getoutput(c)


def preFqs(fastqRoot):
    """
    If the fastq files are well prepared, suitable. 
    """
    fastqs = glob(fastqRoot + "/*.fastq.gz")
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
    for key, fqs in data.items():
        if len(fqs) != 2:
            report = "\t".join([
                str(datetime.now()),
                f"There are no paired fastq files found for %s. Only found %s."
                % (key, ",".join(fqs))
            ])
            print(report)
            del data[key]
    return data


def match(s1, s2, mismatch=2):
    """
    Linker match sequence, allow 2 mismatches.
    """
    s = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            s += 1
    if s > mismatch:
        return False
    else:
        return True


def identifyReads(sample, fin1, fin2, wd, barcodes, linker):
    """
    Identify proper reads. 
    """
    fout1 = wd + "/" + sample + "_R1.fastq.gz"
    fout2 = wd + "/" + sample + "_R2.fastq.gz"
    t = 0
    withBarcode = 0
    withLinker = 0
    i = 0
    with gzip.open(fout1, "wt") as fo1, gzip.open(fout2, "wt") as fo2:
        with gzip.open(fin1, "rt") as f1, gzip.open(fin2, "rt") as f2:
            for r1, r2 in zip(FastqGeneralIterator(f1),
                              FastqGeneralIterator(f2)):
                t += 1
                r1, r2 = list(r1), list(r2)
                #barcode in r2 start position
                b = r2[1][:8]
                if b not in barcodes:
                    continue
                withBarcode += 1
                #check linker sequences
                tmpa = r2[1][8:8 + len(linker)]
                tmpb = r2[1][9:9 + len(linker)]
                tmpc = r2[1][10:10 + len(linker)]
                tmpd = r2[1][11:11 + len(linker)]
                np = -1
                if match(tmpa, linker):
                    np = 8
                if match(tmpb, linker):
                    np = 9
                if match(tmpc, linker):
                    np = 10
                if match(tmpd, linker):
                    np = 11
                #no find of linker
                if np == -1:
                    continue
                #get varialbel C regions
                p = np + len(linker)
                for j in range(np + len(linker), len(r2[1])):
                    if r2[1][j] != "C":
                        p = j
                        break
                withLinker += 1
                #read number id,i5_id, and barcode cell id, to make sure unique id
                rid = "_".join([str(i), sample, b])
                i += 1
                fo1.write("@%s\n%s\n+\n%s\n" % (rid, r1[1], r1[2]))
                fo2.write("@%s\n%s\n+\n%s\n" % (rid, r2[1][p:], r2[2][p:]))
    return fout1, fout2, t, withBarcode, withLinker


def doMapping(sample, fqs, ref, outdir, cpus=25):
    """
    Mapping settings for iscDNase-seq data.
    """
    od = outdir
    sam = od + "/" + sample + ".sam"
    bam = od + "/" + sample + ".bam"
    doBowtie = "bowtie2 --no-mixed --no-discordant --no-unal -p {cpus} -q --local --very-sensitive -x {ref} -1 {fq1} -2 {fq2} -S {sam} ".format(
        cpus=cpus, ref=ref, fq1=fqs[0], fq2=fqs[1], sam=sam)
    report = "\t".join([str(datetime.now()), doBowtie])
    print(report)
    stat, output = subprocess.getstatusoutput(doBowtie)
    #trim with "Warning"
    output = output.split("\n")
    output = [t for t in output if not t.startswith("Warning")]
    output = "\n".join(output)
    report = "\t".join([
        str(datetime.now()), "FLAG_A:" + sample + "\n" + output + "\nFLAG_A\n"
    ])
    print(report)
    with open("iscPre.log", "a") as fo:
        fo.write(report)
    return sample, sam


def sam2bamBedpe(sample, sam, mapq=10):
    """
    SAM to BAM and bedpe file.
    """
    n = os.path.splitext(sam)[0]
    bam = n + ".bam"
    bedpeAll = n + "_all.bedpe"
    #sam to bam, filtering mapq
    samview = "samtools view -b -F 4 -@ 2 -q {mapq} -o {bam} {sam}".format(
        mapq=mapq, bam=bam, sam=sam)
    #sort by read name
    samsort = "samtools sort -n -@ 2 {bam} -T {pre} -o {bam}".format(
        bam=bam, pre=bam.replace(".bam", ""))
    rmsam = "rm %s" % (sam)
    cmds = [samview, samsort, rmsam]
    runCmds(cmds, p=False)
    bam2bedpe = "bamToBed -bedpe -i {bam} > {bedpe}".format(bam=bam,
                                                            bedpe=bedpeAll)
    #report = "\t".join([str(datetime.now()), bam2bedpe])
    #print(report)
    stat, output = subprocess.getstatusoutput(bam2bedpe)
    cmd = "gzip %s" % (bedpeAll)
    runCmds([cmd], p=False)
    return bedpeAll + ".gz"


def sParseBowtie(lines):
    """
    Parse Bowtie2 log file, to obtain mapping stastics.
    """
    d, s = None, None
    lines = lines.split("\n")
    s = lines[0]
    totalReads = int(lines[1].split(";")[0].split()[0])
    mapRatio = float(lines[-2].split("%")[0])
    d = {"TotalRawReads": totalReads, "MappingRatio(%s)": mapRatio}
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


def getStatUniqueTsv(f, fout, blacklist):
    """
    Get unique tsv file and perform statistical analysis. 
    """
    redus = set()
    tot = 0
    unique = 0
    sizes = []
    bs = set()
    #get unique tsv file
    #with gzip.open(fout, "wt") as fo:
    with open(fout + "_tmp", "w") as fo:
        for i, line in enumerate(gzip.open(f, "rt")):
            line = line.split("\n")[0].split("\t")
            if len(line) < 6:
                continue
            tot += 1
            rid = line[6]
            #record barcode number
            b = "_".join(rid.split("_")[1:])
            bs.add(b)
            if line[0] != line[3]:
                continue
            chrom = line[0]
            start = min([int(line[1]), int(line[4])])
            end = max([int(line[2]), int(line[5])])
            #remove redudant PETs
            r = hash((chrom, start, end, b))
            if r in redus:
                continue
            else:
                redus.add(r)
            unique += 1
            line = [chrom, str(start), str(end), b]
            fo.write("\t".join(line) + "\n")
            sizes.append(end - start)
    c1 = f"intersectBed -a {fout}_tmp -b {blacklist} -v | gzip > {fout}"
    c2 = f"rm {fout}_tmp"
    runCmds([c1, c2], p=False)
    finalN = 0
    for line in gzip.open(fout, "rt"):
        finalN += 1
    return tot, unique, finalN, np.median(sizes), len(bs)


@click.command()
@click.option(
    "-d",
    required=True,
    type=str,
    help=
    "The directory for raw .fastq.gz files, for example ../1.fastq/. Read 1 ends with _R1.fastq.gz and read 2 ends with _R2.fastq.gz. "
)
@click.option(
    "-o",
    required=True,
    help="Output prefix",
    type=str,
)
@click.option(
    "-barcode",
    required=True,
    type=str,
    help="Barcode information for marking single-cell, one line one barcode.")
@click.option(
    "-blacklist",
    required=True,
    type=str,
    help=
    "ENCODE blacklist file in BED format to filtering un-reliable reads. Without filtering, aggregated signals will be biased."
)
@click.option(
    "-ref",
    required=True,
    type=str,
    help=
    "Bowtie2 reference index prefix, such as ./ref/hg38, generated from bowtie2-build hg38.fa hg38."
)
@click.option("-linker",
              required=False,
              type=str,
              default="AGAACCATGTCGTCAGTGT",
              help="Linker sequence after barcode in R2.")
@click.option("-p",
              default=4,
              type=int,
              help="Number of CPUs to finish the job, default is set to 4.")
def run(d, o, barcode, blacklist, ref, linker, p=4, mapq=10):
    """
    Pre-process raw reads from indexing split and pool single-cells data to genomic regions.

    Example: 

    python iscPre.py -d ../test/test -o test -barcode barcodes.txt -blacklist encode_mm10_blacklist.bed -ref bowtie2_mm10 -p 15 
    """
    #step 0, prepare everything
    for t in ["bowtie2", "samtools", "bamToBed", "intersectBed"]:
        if not isTool(t):
            report = "\t".join([
                str(datetime.now()),
                f"{t} not exits! Please install through conda."
            ])
            print(report)
            return
    if not os.path.exists(d):
        report = "\t".join(
            [str(datetime.now()), f"Input {d} not exists! Return."])
        print(report)
        return
    if len(glob(ref + "*.bt2")) == 0:
        report = "\t".join([
            str(datetime.now()),
            f"Bowtie2 reference not exists for prefix of {ref}! Return."
        ])
        print(report)
        return
    if not os.path.exists(o):
        os.makedirs(o, exist_ok=True)
    else:
        fs = glob(os.path.join(o, "*"))
        if len(fs) > 0:
            report = "\t".join([
                str(datetime.now()),
                f"Target output directory {o} is not empty, may over-write some files."
            ])
            print(report)
    if not os.path.isfile(barcode):
        report = "\t".join(
            [str(datetime.now()), f"Barcode file {barcode} not exists."])
        print(report)
    if not os.path.isfile(blacklist):
        report = "\t".join(
            [str(datetime.now()), f"Barcode file {blacklist} not exists."])
        print(report)

    #working directory for tmp files
    wd = str(random.random())
    os.mkdir(wd)

    #step 1, trim read2 and reid all reads
    report = "\t".join([
        str(datetime.now()),
        "Step1: Trim linkers and rename read with barcodes."
    ])
    print(report)

    #read in barcodes
    barcodes = open(barcode).read().split("\n")[:-1]
    barcodes = [b.strip() for b in barcodes]
    barcodes = set(barcodes)
    #seperate fastq files
    fqs = preFqs(d)
    #identify proper reads
    fout1 = o + "/" + o + "_R1.fastq.gz"
    fout2 = o + "/" + o + "_R2.fastq.gz"
    ds = Parallel(n_jobs=p)(
        delayed(identifyReads)(sample, fqs[0], fqs[1], wd, barcodes, linker)
        for sample, fqs in tqdm(fqs.items()))

    report = "\t".join(
        [str(datetime.now()), "Step 2: Merge all seperate fastq files."])
    print(report)

    fq1s = []
    totalRaw, rawWithBarcode, rawWithLinker = 0, 0, 0
    for d in tqdm(ds):
        cmd1 = "cat %s >> %s" % (d[0], fout1)
        cmd2 = "cat %s >> %s" % (d[1], fout2)
        runCmds([cmd1, cmd2], p=False)
        totalRaw += d[2]
        rawWithBarcode += d[3]
        rawWithLinker += d[4]
    shutil.rmtree(wd)

    #step 2, mapping
    report = "\t".join(
        [str(datetime.now()), "Step 3: Map processed reads to genome."])
    print(report)
    sample, sam = doMapping(o, [fout1, fout2], ref, o, cpus=p)

    #step 3, convert to bam and bedpe files
    #sam to bam and bedpe
    report = "\t".join([str(datetime.now()), "Step 4: File type conversion."])
    print(report)
    bedpe = sam2bamBedpe(sample, sam, mapq)

    #step 4, get unique tsv
    report = "\t".join([str(datetime.now()), "Step 5: Analyze final reads."])
    print(report)
    tsvUni = o + "/" + o + "_unique.tsv.gz"
    totalHighQuality, uniqueHighQuality, uniqueHighQaulityRemoveBlacklist, medianSize, barcodeNumber = getStatUniqueTsv(
        bedpe, tsvUni, blacklist)

    #step 5, summarize everything
    mata = parseBowtielog()
    if totalHighQuality > 0:
        redundancy = 1 - float(uniqueHighQuality) / totalHighQuality
    else:
        redundancy = 0
    summary = {
        "1_totalRawReads": totalRaw,
        "2_rawReadsWithBarcode": rawWithBarcode,
        "3_rawReadsWithBarcodeLinker": rawWithLinker,
        "4_mappingRatio": mata.loc[o, "MappingRatio(%s)"] / 100,
        "5_mappedReads(mapq>=%s)" % mapq: totalHighQuality,
        "6_mappedHighQualityUniqueReads": uniqueHighQuality,
        "7_redundancy": redundancy,
        "8_finalUniqueReadsRemoveBlacklist": uniqueHighQaulityRemoveBlacklist,
        "9_medianFragmentSize": medianSize,
        "10_barcodeNumber": barcodeNumber,
        "11_yield": float(uniqueHighQaulityRemoveBlacklist)/totalRaw,
    }
    summary = pd.Series(summary)
    fo = o + "_iscPreFq_summary.txt"
    summary.to_csv(fo, sep="\t", header=None)

    #finish
    report = "\t".join([
        str(datetime.now()),
        f"Processing of {o} finished. Summary written to {fo}."
    ])
    print(report)


if __name__ == '__main__':
    run()
