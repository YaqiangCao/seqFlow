#!/usr/bin/env python2.7
#--coding:utf-8 --
"""
peaksOverlapEnrich.py
The aim is to detect  the significant overlaps between target sets and annotation sets based on fisher exact test. However, fisher exact test needs the background, so except the file for test, all others are used as background, which can be see to choose the relative significant annotation. 
2019-11-25: updated by scipy.stat for significance test. 
"""

__author__ = "CAO Yaqiang"
__date__ = "2016-03-22"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import os, copy
from glob import glob
from collections import Counter

#3rd library
import pandas as pd
import numpy as np
from tqdm import tqdm
from scipy.stats import fisher_exact
from scipy.stats import chi2_contingency as chi2
from joblib import Parallel, delayed

#my setting
from utils import *


def getCov(bedf):
    print("building coverage for %s" % bedf)
    cov = {}
    number = 0
    bps = 0
    for line in open(bedf):
        number += 1
        line = line.split("\n")[0].split("\t")
        c = line[0]
        s = int(line[1])
        e = int(line[2])
        bps += e - s
        if c not in cov:
            cov[c] = set()
        for i in range(s, e + 1):
            cov[c].add(i)
    return number, bps, cov


def getBgCov(fs):
    print("building background coverage")
    cov = {}
    for f in tqdm(fs):
        for line in open(f):
            line = line.split("\n")[0].split("\t")
            c = line[0]
            s = int(line[1])
            e = int(line[2])
            if c not in cov:
                cov[c] = set()
            for i in range(s, e + 1):
                cov[c].add(i)
    number = 0
    bps = 0
    #get the merged sets of regions
    for c in cov.keys():
        s = list(cov[c])
        s.sort()
        i = 0
        while i < len(s):
            j = i+1
            while j < len(s):
                if s[j] - s[j-1] > 1:
                    break
                else:
                    j += 1
            j = j -1
            number += 1
            bps += s[j] - s[i]
            i = j+1
    return number, bps, cov


def getStat(number, bps, cov, rs):
    overlaps = 0
    overlapBps = 0
    bpsb = 0
    for r in rs:
        flag = False
        bpsb += r[2] - r[1]
        if r[0] not in cov:
            continue
        for i in range(r[1], r[2]):
            if i in cov[r[0]]:
                flag = True
                overlapBps += 1
        if flag == True:
            overlaps += 1
    jinumber = float(overlaps) / (number + len(rs) - overlaps)
    jibp = float(overlapBps) / (bps + bpsb - overlapBps)
    #print(number, len(rs), overlaps, overlaps / len(rs), bps, bpsb, overlapBps / bpsb, jinumber, jibp)    
    return overlaps, overlapBps, jinumber, jibp, bpsb


def getEnrich(targetf, anofs):
    rs = []
    for line in open(targetf):
        line = line.split("\n")[0].split("\t")
        if len(line) < 3:
            continue
        line[1] = int(line[1])
        line[2] = int(line[2])
        rs.append(line)
    ds = {}
    for f in anofs:
        name = f.split("/")[-1].split(".")[0]
        nfs = [nf for nf in anofs if nf != f]
        number, bps, cov = getCov(f)
        numberBg, bpsBg, covBg = getBgCov(nfs)
        overlaps, overlapBps, jinumber,jibp, bpsb = getStat( number,bps,cov,rs  )
        overlapsBg, overlapBpsBg, jinumberbg, jibpbg, bpsb = getStat( numberBg,bpsBg,covBg,rs  )
        oddsratio, pvalue = fisher_exact([[ overlaps, overlapsBg], [number - overlaps, numberBg-overlapsBg]])  
        oddsratio_bp, pvalue_bp = fisher_exact([[ overlapBps, overlapBpsBg], [bps - overlapBps, bpsBg-overlapBpsBg]])  
        chisq, p = chi2([
                            [ overlaps, number-overlaps,],
                            [ overlapsBg, numberBg-overlapsBg]
                        ])[:2]
        chisq, pbp = chi2([
                            [ overlapBps, bps-overlapBps,],
                            [ overlapBpsBg, bpsBg-overlapBpsBg]
                        ])[:2]
        ds[name] = {
            "regions": len(rs),
            "tRegions": number,
            "bgRegions": numberBg,
            "regionBps": bpsb,
            "tRegionBps": bps,
            "bgRegionBps": bpsBg,
            "tOverlaps": overlaps,
            "bgOverlaps":overlapsBg,
            "tJInumber": jinumber,
            "tJIbp": jibp,
            "bgJInumber": jinumberbg ,
            "bgJIbp": jibpbg,
            "oddRatioNumber": oddsratio,
            "fisherPvalueNumber":pvalue,
            "oddRatioBp": oddsratio_bp,
            "fisherPvalueBp":pvalue_bp,
            "chiPvalue":p,
            "chiPvalueBp":pbp,
            }
        print(n,ds[name])
    ds = pd.DataFrame(ds)
    ds.to_csv("stat.txt",sep="\t")


def main():
    targetf = "../../2.MergedAcnhors/GM12878_Trac_anchors.bed"
    anofs = glob(
        "/home/caoy7/caoy7/Projects/0.Reference/1.hg38/6.ChIPPeaks/2.remap/GM12878/*.bed"
    )
    anofs.sort()
    getEnrich(targetf, anofs)


if __name__ == "__main__":
    main()
