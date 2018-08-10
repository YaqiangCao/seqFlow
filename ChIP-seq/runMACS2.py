#!/usr/bin/env python2.7
#--coding:utf-8--
"""
repeatom_macs2_peakcalling.py
Basicially using the default parameters.
"""

#general library
import glob, string, argparse, os, os.path, commands, time, logging, sys
from datetime import datetime

#3rd library
import HTSeq
from joblib import Parallel, delayed

#own library
from Biolibs.rel.General.logger import getlogger
from Biolibs.rel.General.countLines import count_lines
date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
logger = getlogger(fn=os.getcwd() + "/" + date.strip() + "_" +
                   os.path.basename(__file__) + ".log")

__author__ = "CAO Yaqiang"
__date__ = "2014-11-13"
__modified__ = "2014-11-14"
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


def peakCalling(cbed, tbed, subpeak=True, root="1.PeakCalling"):
    sample = os.path.splitext(os.path.split(tbed)[1])[0]
    peak_f = root + "/" + sample + "/" + sample + "_peaks.bed"
    if os.path.exists(peak_f):
        print peak_f, "exists!"
        return
    reads_cbed = count_lines(cbed)
    reads_tbed = count_lines(tbed)
    doMACS = "macs2 callpeak -c {cbed} -t {tbed} -n {sample} -g hs --nomodel --extsize 73 --keep-dup 1 --outdir {out}".format(
        cbed=cbed, tbed=tbed, sample=sample, out=root + "/" + sample)
    #in normal situation, input reads should larger than tf specific reads, which is macs default set (common happens in public available data.)
    if reads_tbed > reads_cbed:
        doMACS += " --to-large"
    if subpeak:
        doMACS += " --call-summits"
    call_sys([doMACS])
    logger.info("%s finished PeakCalling!" % sample)


def extendBed_to_bed(f):
    model = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    #collect peaks
    for line in open(f):
        line = line.split("\n")[0].split("\t")
        if len(line) < 3:
            continue
        try:
            iv = HTSeq.GenomicInterval(line[0],
                                       int(line[1]), int(line[2]), "+")
        except:
            continue
        peak_name = line[3]
        if peak_name[-1] in string.letters:
            peak_name = peak_name[:-1]
        model[iv] += peak_name
    #re-orgonize
    peaks = []
    for iv, value in model.steps():
        if value:
            p = [iv.chrom, str(iv.start), str(iv.end), list(value)[0]]
            peaks.append(p)
    return peaks


def prepare_samples(root="../3.Mapping/"):
    bams = glob.glob(root + "*/*.bam")
    samples = {}
    for b in bams:
        n = os.path.splitext(os.path.split(b)[1])[0]
        cell = "_".join(n.split("_")[:-1])
        target = n.split("_")[-1]
        if cell not in samples:
            samples[cell] = {"input": None, "target": []}
        if target == "Input":
            samples[cell]["input"] = b
        else:
            samples[cell]["target"].append(b)
    nsamples = []
    for s in samples:
        for t in samples[s]["target"]:
            nsamples.append((samples[s]["input"], t))
    return nsamples


def main():
    root = "../3.Mapping/"
    samples = prepare_samples(root)
    root = "1.PeakCalling"
    if not os.path.exists(root):
        os.mkdir(root)
    Parallel(n_jobs=11)(delayed(peakCalling)(b[0], b[1]) for b in samples)
    fs = glob.glob(root + "/*/*narrowPeak")
    fs.extend(glob.glob(root + "/*/*broadPeak"))
    for f in fs:
        peaks = extendBed_to_bed(f)
        fn = os.path.splitext(f)[0] + ".bed"
        print fn
        with open(fn, "w") as f2:
            for p in peaks:
                f2.write("\t".join(p) + "\n")


if __name__ == "__main__":
    startTime = datetime.now()
    main()
    elapsed = datetime.now() - startTime
    print "The process is done."
    print "Time used:", elapsed
