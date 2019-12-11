#!/usr/bin/env python2.7
#--coding:utf-8--
"""
cutTracFq.py
2019-12-10:
"""

__author__ = "CAO Yaqiang"
__date__ = "2019-12-10"
__email__ = "caoyaqiang0410@gmail.com"

#sys
import gzip,time,os
from glob import glob
from datetime import datetime

#3rd
from joblib import Parallel, delayed
from Bio.SeqIO.QualityIO import FastqGeneralIterator

#global settings
TARGET="CTGTCTCTT"


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


def filterSeq(fq1,fq2,pre,lencut=15): 
    """
    Filtering target adapter/linker seq in sequencing data.
    """
    fout1 = gzip.open(pre + "_1.fastq.gz","wt")
    fout2 = gzip.open(pre + "_2.fastq.gz","wt")
    r1a = 0 
    r2a = 0
    short= 0
    total = 0
    #with open(fq1) as f1, open(fq2) as f2:
    with gzip.open(fq1,"rt") as f1, gzip.open(fq2,"rt") as f2:
        for r1, r2 in zip(FastqGeneralIterator(f1), FastqGeneralIterator(f2)):
            total += 1
            r1, r2 = list(r1), list(r2)
            if TARGET in r1[1]:
                p = r1[1].find(TARGET) + len(TARGET)
                r1[1] = r1[1][p+1:]
                r1[2] = r1[2][p+1:]
                r1a +=1 
            else:
                r1[1] = r1[1][5:]
                r1[2] = r1[2][5:]
            if TARGET in r2[1]:
                p = r2[1].find(TARGET) + len(TARGET)
                r2[1] = r2[1][p+1:]
                r2[2] = r2[2][p+1:]
                r2a += 1
            else:
                r2[1] = r2[1][5:]
                r2[2] = r2[2][5:]
            if len(r1[1])< lencut or len(r2[1]) < lencut:
                short += 1
                continue
            fout1.write("@%s\n%s\n+\n%s\n" % (r1[0], r1[1], r1[2]))
            fout2.write("@%s\n%s\n+\n%s\n" % (r2[0], r2[1], r2[2]))
    print("total %s pairs;\nfor read1, %.3f has adapter;\nfor read2, %.3f has adapter;\n%.3f pairs too short (any of them < 15bp)"%(total, float(r1a)/total,float(r2a)/total, float(short)/total ))
    fout1.close()
    fout2.close()
    return sample,total, float(r1a)/total,float(r2a)/total,float(short)/total


def main():
    data = prepare_fastq("../../1.fastq/")
    jobs = min(10,len(data))
    ds = Parallel(n_jobs=10)(delayed(filterSeq)(fqs[0],fqs[1],sample )
                       for sample, fqs in data.items())
    data = {}
    for d in ds:
        data[ d[0] ] = {
            "total":d[1],
            "r1_adapter_ratio":d[2],
            "r2_adapter_ratio":d[3],
            "filtered_short_ratio":d[4],
            }
    data = pd.DataFrame(data).T
    data.to_csv("filterAdapterStat.txt",sep="\t")


if __name__ == '__main__':
    start_time = datetime.now()
    date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
    main()
    elapsed = datetime.now() - start_time
    print("The process is done")
    print("Time used:", elapsed)
