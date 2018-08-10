#!/usr/bin/env python2.7
#--coding:utf-8 --
"""
peak2geneTSS.py
2015-01-11: Small bug for shift fixed.
"""

__author__ = "CAO Yaqiang"
__date__ = "2014-09-29"
__modified__ = "2014-11-25"
__version__ = "0.1"
__email__ = "caoyaqiang0410@gmail.com"

#general library
import sys, glob, copy
try:
    import cPickle as pickle
except:
    import pickle

#3rd library
import HTSeq
import pandas as pd
import numpy as np
from joblib import Parallel, delayed


def get_TSS_TES(gtf):
    """
    """
    regions = {}
    for g in HTSeq.GFF_Reader(gtf):
        if "gene_name" not in g.attr:
            continue
        gid = g.attr["gene_name"]
        if ":" in gid:
            gid = gid.split(":")[1]
        if gid not in regions:
            regions[gid] = g.iv
        else:
            if g.iv.start < regions[gid].start:
                regions[gid].start = g.iv.start
            if g.iv.end > regions[gid].end:
                regions[gid].end = g.iv.end
    return regions


def shift_TSS_TES(regions, shift=[2000, 2000]):
    """
    Strand specific extension., Specific for TSS.
    """
    for key, iv in regions.items():
        #only use the TSS
        iv.end = iv.start
        if iv.strand == "+":
            if iv.start - shift[0] < 0:
                iv.start = 0
            else:
                iv.start = iv.start - shift[0]
            iv.end = iv.end + shift[1]
        else:
            if iv.end - shift[1] < 0:
                iv.end = 0
            else:
                iv.end = iv.start - shift[1]
            iv.start = iv.start + shift[0]
        regions[key] = iv
    return regions


def get_gene_model(gtf=None, shift=[2000, 2000]):
    bgmodel = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    regions = get_TSS_TES(gtf)
    regions = shift_TSS_TES(regions, shift=shift)
    for key, iv in regions.items():
        if iv.start > iv.end:
            # here the model is nonstranded 
            iv.start, iv.end = iv.end, iv.start
        bgmodel[iv] += key
    return bgmodel


def get_gene(bgmodel=None, bed=None):
    gs = set()
    for item in HTSeq.BED_Reader(bed):
        iv = item.iv
        #just in case the bed contains not included chrom in background model
        try:
            for ivb, valueb in bgmodel[iv].steps():
                gs.update(valueb)
        except:
            continue
    return gs


def get_annotations(beds, fnOut, gtf):
    bgmodel = get_gene_model(gtf=gtf, shift=[2000, 2000])
    data = {}
    for bed in beds:
        tf = bed.split("/")[-1].split(".")[0]
        gs = get_gene(bgmodel=bgmodel, bed=bed)
        if len(gs) == 0:
            continue
        if tf not in data:
            data[tf] = gs
        else:
            data[tf].update(gs)
    with open(fnOut, "w") as f:
        for key, gs in data.items():
            gs = list(gs)
            gs = map(str, gs)
            gs = [g for g in gs if g != "NaN" and g != "nan"]
            line = key + "\t" + ",".join(gs) + "\n"
            f.write(line)


def main():
    gtf = "genecode.v21.appris.pcRNA.gtf"
    beds = glob.glob("../1.Classfiled_Merged/3.conservedFactors/*.bed")
    fnOut = "ENCODE_mRNA_TSS_conserved.gmt"
    get_annotations(beds, fnOut, gtf)
    beds = glob.glob("../1.Classfiled_Merged/1.mergeReps/*.bed")
    fnOut = "ENCODE_mRNA_TSS_all.gmt"
    get_annotations(beds, fnOut, gtf)
    beds = glob.glob("../1.Classfiled_Merged/2.mergeFactors/*.bed")
    fnOut = "ENCODE_mRNA_TSS_factors.gmt"
    get_annotations(beds, fnOut, gtf)
    gtf = "genecode.v21.gtf"
    beds = glob.glob("../1.Classfiled_Merged/3.conservedFactors/*.bed")
    fnOut = "ENCODE_all_TSS_all_conserved.gmt"
    get_annotations(beds, fnOut, gtf)
    beds = glob.glob("../1.Classfiled_Merged/1.mergeReps/*.bed")
    fnOut = "ENCODE_all_TSS_all.gmt"
    get_annotations(beds, fnOut, gtf)
    beds = glob.glob("../1.Classfiled_Merged/2.mergeFactors/*.bed")
    fnOut = "ENCODE_all_TSS_factors.gmt"
    get_annotations(beds, fnOut, gtf)


if __name__ == "__main__":
    main()
