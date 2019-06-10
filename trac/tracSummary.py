#!/usr/bin/env python2.7
#--coding:utf-8--
"""
2019-05-21
"""

__author__ = "CAO Yaqiang"
__date__ = ""
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import os, time
from glob import glob
from datetime import datetime

#3rd library
import numpy as np
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed


def getReport(mapf, bedpef):
    mapd = pd.read_table(mapf, index_col=0, sep="\t")
    bd = pd.read_table(bedpef, index_col=0, sep="\t")
    data = mapd.join(bd)
    print(list(data.columns))
    data["redundancy(%s)"] = (
        1.0 - data["NoneRedudantPETs"] /
        (data["TotalReads"] * data["MappingRatio(%s)"] / 100)) * 100
    data["yield(%s)"] = data["cisPETs"] / data["TotalReads"] * 100
    data = data[[
        "TotalReads", "ConcordantlyUniqueMapReads",
        "DisconcordantlyUniqueMapReads", "MappingRatio(%s)",
        "NoneRedudantPETs", "redundancy(%s)", "cisPETs", "yield(%s)",
        "distance<=150bp", "150<=distance<=1000bp", "1kb<distance<=10k",
        "distance>20kb"
    ]]
    data.to_csv("qualityReport.txt", sep="\t")


def main():
    getReport("../3.mapping/MappingStat.txt",
              "../4.uniqueNonRedudantBEDPE/3.evaBedpe/bedpe_quality_stat.txt")


if __name__ == '__main__':
    start_time = datetime.now()
    main()
    elapsed = datetime.now() - start_time
    print "The process is done"
    print "Time used:", elapsed
