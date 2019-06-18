#!/usr/bin/env python2.7
#--coding:utf-8--
"""
callVariants.py
2018-08-08
"""

__author__ = "CAO Yaqiang"
__date__ = "2018-07-04"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import os, time
from glob import glob
from datetime import datetime

#3rd library
import numpy as np
import pandas as pd
import seaborn as sns
from joblib import Parallel, delayed
#my own
from Biolibs.rel.General.logger import getlogger

#global settings
#logger
date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
logger = getlogger(fn=os.getcwd() + "/" + date.strip() + "_" +
                   os.path.basename(__file__) + ".log")

#global
#v4 not compatitable for RNA-seq analysis pipeline
#gatk="/picb/molsysbio2/caoyaqiang/2.WWX_RNA-seq/3.BigSamples/7.mutations/gatk-4.0.5.2/gatk"
#roll back to v3
gatk = "/picb/molsysbio2/caoyaqiang/2.WWX_RNA-seq/3.BigSamples/7.mutations/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar"
picard = "/picb/molsysbio2/caoyaqiang/2.WWX_RNA-seq/3.BigSamples/7.mutations/picard.jar"
snpsift = "/picb/molsysbio2/caoyaqiang/2.WWX_RNA-seq/0.Tools/snpEff/SnpSift.jar"
FA = "/home/caoyaqiang/cyq_m2/2.WWX_RNA-seq/1.Reference/1.FAs/hg38.fa"
dbSNP = "/picb/molsysbio2/caoyaqiang/2.WWX_RNA-seq/1.Reference/5.SNPs/dbSNP_all_v151.vcf.gz"


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


def prepareBam(f):
    root = "1.PreReads/"
    n = f.split("/")[-1].split("_Aligned")[0]
    fo = root + n + ".bam"
    mark = root + n + ".mark"
    ta = root + n + ".tmpa.bam"
    tb = root + n + ".tmpb.bam"
    if os.path.isfile(fo) or os.path.isfile(mark):
        logger.info("%s has been processed" % fo)
        return
    with open(mark, "w") as ft:
        ft.write("\n")
        ft.close()
    cmd1 = "java -Djava.io.tmpdir=./tmp -jar {picard} AddOrReplaceReadGroups I={input} O={output} SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample".format(
        picard=picard, input=f, output=ta)
    cmd2 = "java -Djava.io.tmpdir=./tmp -jar {picard} MarkDuplicates I={input} O={output} CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics".format(
        picard=picard, input=ta, output=tb)
    cmd3 = "java -Djava.io.tmpdir=./tmp -jar {gatk} -T SplitNCigarReads -R {fa} -I {input} -o {output} -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS".format(
        gatk=gatk, input=tb, output=fo, fa=FA)
    cmd4 = "rm %s %s %s" % (ta, tb, mark)
    call_sys([cmd1, cmd2, cmd3, cmd4])
    return


def callMutation(f):
    root = "2.VariantsCalling/"
    n = f.split("/")[-1].split(".bam")[0]
    fo = root + n + ".vcf"
    if os.path.exists(fo):
        return
    cmd1 = "java -jar {gatk} -T HaplotypeCaller -R {FA} -I {input} -o {output} -dontUseSoftClippedBases -stand_call_conf 20.0".format(
        gatk=gatk, FA=FA, input=f, output=fo)
    call_sys([cmd1])


def filterMutation(f):
    if not os.path.isfile(f + ".idx"):
        return
    root = "3.VariantsFiltering/"
    fo = root + f.split("/")[-1]
    if os.path.exists(fo):
        return
    cmd = 'java -jar {gatk} -T VariantFiltration -R {FA} -V {input} -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o {output}'.format(
        gatk=gatk, FA=FA, input=f, output=fo)
    call_sys([cmd])


def bgzipTabix(f):
    if not f.endswith(".vcf"):
        return
    cmd1 = "bgzip %s" % f
    cmd2 = "bcftools index %s.gz" % f
    cmd3 = "tabix %s.gz" % f
    call_sys([cmd1, cmd2, cmd3])


def anoVcf(f):
    root = "4.VariantsAnnotation/"
    n = f.split("/")[-1].split(".vcf.gz")[0]
    ft = root + n + ".temp"
    fo = root + n + ".vcf"
    if os.path.isfile(ft) or os.path.isfile(fo) or os.path.isfile(fo + ".gz"):
        return
    chroms = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM"
    cmd1 = 'bcftools filter -O v -r {chroms} {fin} | grep -v "contig" | grep -v "ERCC" > {fo}'.format(
        chroms=chroms, fin=f, fo=ft)
    cmd2 = 'java -jar {snpsift} annotate {dbSNP} {fin} > {fo}'.format(
        snpsift=snpsift, dbSNP=dbSNP, fin=ft, fo=fo)
    cmd3 = "rm %s" % ft
    call_sys([cmd1, cmd2, cmd3])
    bgzipTabix(fo)


def main():
    #fs = glob("../2.Mapping/*/*.bam")
    #Parallel( n_jobs=8 )( delayed( prepareBam )( bam ) for bam in fs)
    #fs = glob("1.PreReads/*.bam")
    #Parallel( n_jobs=5 )( delayed( callMutation )( bam ) for bam in fs)
    #fs = glob("2.VariantsCalling/*.vcf")
    #Parallel( n_jobs=10 )( delayed( filterMutation )( vcf ) for vcf in fs)
    #fs = glob("2.VariantsCalling/*.vcf")
    #fs.extend(glob("3.VariantsFiltering/*.vcf"))
    #Parallel( n_jobs=20 )( delayed( bgzipTabix )( vcf ) for vcf in fs)
    fs = glob("2.VariantsCalling/*.vcf.gz")
    Parallel(n_jobs=10)(delayed(anoVcf)(vcf) for vcf in fs)


main()
