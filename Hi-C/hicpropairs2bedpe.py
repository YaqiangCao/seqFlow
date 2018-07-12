#!/usr/bin/env pyhthon2.7
#--coding:utf-8--
"""
# Author: Ai Daosheng
# Contact: remyai@hotmail.com
# Created Time : Mon 25 Sep 2017 04:10:24 PM CST
# Last Modified :
# File Name: hicpropairs2bedpe.py
# Description: convert pets format from HiC-Pro to bedpe format

"""

import os,gzip
from glob import glob 
from multiprocessing import Pool

cs = ["chr%s"%i for i in xrange(1,23)]
cs.extend(["chrX","chrY","chrM"])

def pairs2bedpe(f_hicpro):
    print f_hicpro
    with gzip.open(f_hicpro.replace("allValidPairs","bedpe"),'w') as f_bedpe:
    #with open(f_hicpro.replace("allValidPairs","bedpe"),'w') as f_bedpe:
    #with open(f_hicpro.replace("validPairs","bedpe"),'w') as f_bedpe:
        with gzip.open(f_hicpro) as f_pair:
        #with open(f_hicpro) as f_pair:
            for line in f_pair:
                line = line.strip().split('\t')
                if line[3] == '+':
                    petA = [line[1],line[2],int(line[2])+int(line[-2])]
                else:
                    petA = [line[1],int(line[2])-int(line[-2]),line[2]]
                if line[6] == '+':
                    petB = [line[4],line[5],int(line[5])+int(line[-1])]
                else:
                    petB = [line[4],int(line[5])-int(line[-1]),line[5]]
                newline = [petA[0],petA[1],petA[2],petB[0],petB[1],petB[2],line[0],'.',line[3],line[6]]
                if newline[0] != newline[3]:
                    continue
                if newline[0] not in cs:
                    continue
                f_bedpe.write('\t'.join(map(str,newline))+'\n')
#map(pairs2bedpe,glob("*allValidPairs"))
map(pairs2bedpe,glob("hicpro/hic_results/data/*/*allValidPairs"))
