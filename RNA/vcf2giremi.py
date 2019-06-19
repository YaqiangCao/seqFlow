#!/usr/bin/env python2.7
#--coding:utf-8--
"""
vcf2gerimiL.py
2018-08-08
"""

__author__ = "CAO Yaqiang"
__date__ = "2018-08-08"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import os, gzip
from glob import glob

#3rd library
import HTSeq
from joblib import Parallel, delayed


def getGenes():
    bed = "/picb/molsysbio2/caoyaqiang/2.WWX_RNA-seq/1.Reference/2.Genecode/2.Parsed/genecode.v27.bed"
    genes = {}
    for line in open(bed):
        line = line.split("\n")[0].split("\t")
        attrs = {}
        for t in line[9].split(";"):
            t = t.split(" ")
            if len(t) == 1 or '"' not in t[1]:
                continue
            t[1] = t[1].split('"')[1]
            attrs[t[0]] = t[1]
        gid = attrs["gene_id"] + "|" + attrs["gene_name"]
        s = int(line[1])
        e = int(line[2])
        if gid not in genes:
            genes[gid] = {
                "chr": line[0],
                "start": s,
                "end": e,
                "strand": line[5]
            }
        else:
            if s < genes[gid]["start"]:
                genes[gid]["start"] = s
            if e > genes[gid]["end"]:
                genes[gid]["end"] = e
    return genes


def getGeneModel():
    genes = getGenes()
    model = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    for g in genes.keys():
        iv = HTSeq.GenomicInterval(genes[g]["chr"],
                                   genes[g]["start"],
                                   genes[g]["end"],
                                   strand=genes[g]["strand"])
        model[iv] += g
    return model, genes


def convert(vcfIn):
    fout = vcfIn.split("/")[-1].split(".vcf.gz")[0] + ".gerimi"
    if os.path.isfile(fout):
        return
    model, genes = getGeneModel()
    with open(fout, "w") as fo:
        for line in gzip.open(vcfIn):
            if line.startswith("#"):
                continue
            line = line.split("\n")[0].split("\t")
            p = HTSeq.GenomicPosition(line[0], int(line[1]))
            r = list(model[p])
            if len(r) == 0:
                g = "Inte"
                strand = "#"
            else:
                g = r[0]
                strand = genes[g]["strand"]
            if line[2].startswith("rs"):
                flag = 1
            else:
                flag = 0
            nline = [line[0], line[1], int(line[1]) + 1, g, flag, strand]
            fo.write("\t".join(map(str, nline)) + "\n")


def main():
    fs = glob("../../7.mutations/4.VariantsAnnotation/*.vcf.gz")
    Parallel(n_jobs=10)(delayed(convert)(vcf) for vcf in fs)


main()
