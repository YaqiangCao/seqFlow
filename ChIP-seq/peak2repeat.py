#!/usr/bin/env python2.7
#--coding:utf-8--
"""
peaks2repeats.py
"""

#general library
import glob, string, os, time, logging, sys
from datetime import datetime

#3rd library
import HTSeq
import pandas as pd
from joblib import Parallel, delayed

REPEATS_LOCUS = "/picb/molsysbio/usr/caoyaqiang/1.Projects/17.Alu/1.GenomeReference/1.hg38/6.Repeats_Annotation/1.GenomicLocation/hg38_rep_locus.txt"

__author__ = "CAO Yaqiang"
__date__ = "2015-01-28"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"


def read_rep_locus(shift=[0, 0]):
    rep_model = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    for line in open(REPEATS_LOCUS):
        line = line.split("\n")[0].split("\t")
        if len(line) != 6:
            continue
        locus = [line[0], line[1], line[2], line[3]]
        rep = "|".join(line)
        iv = HTSeq.GenomicInterval(
            locus[0], int(locus[1]), int(locus[2]), strand=locus[3])
        if iv.start - shift[0] < 0:
            iv.start = 0
        else:
            iv.start = iv.start - shift[0]
        iv.end = iv.end + shift[1]
        rep_model[iv] += rep
    return rep_model


def get_overlap_rep(rep_model, peak_f):
    reps = set()
    c_t = os.path.split(peak_f)[1].split("_")[0]
    for g in HTSeq.BED_Reader(peak_f):
        iv = g.iv
        try:
            for niv, value in list(rep_model[iv].steps()):
                reps.update(value)
        except:
            continue
    return c_t, reps


def get_binary_matrix(rep_model, pre):
    peaks = glob.glob("../1.Peaks/*.bed")
    data = {}
    for p in peaks:
        print p
        c_t, reps = get_overlap_rep(rep_model, p)
        print c_t, " read finished!"
        data[c_t] = {}
        for r in reps:
            data[c_t][r] = 1
    data = pd.DataFrame(data)
    data = data.fillna(0)
    nis = []
    for i in data.index:
        if data.loc[i, :].sum() == 0:
            nis.append(i)
    data = data.drop(nis)
    fn = pre + ".txt"
    data.to_csv(fn, sep="\t", index_label="rep")
    print "finished!"


def main():
    for d in [0, 1, 2, 5, 10]:
        pre = "DHS_Shift_%sk" % d
        rep_model = read_rep_locus(shift=[1000 * d, 1000 * d])
        get_binary_matrix(rep_model, pre)


if __name__ == "__main__":
    startTime = datetime.now()
    main()
    elapsed = datetime.now() - startTime
    print "The process is done."
    print "Time used:", elapsed
