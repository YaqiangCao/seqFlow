#!/usr/bin/env python2.7
#--coding:utf-8--
"""
starMapping.py
"""

__author__ = "CAO Yaqiang"
__date__ = "2014-10-15"
__modified__ = "2015-03-02"
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import glob, os, time, commands, shutil
from datetime import datetime

#3rd library
import matplotlib as mpl
mpl.use("pdf")
import seaborn as sns
import pylab
sns.set_style("white")
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["figure.figsize"] = (4, 2.75)
mpl.rcParams["figure.dpi"] = 100
mpl.rcParams["savefig.transparent"] = True
mpl.rcParams["savefig.bbox"] = "tight"
mpl.rcParams["font.size"] = 10.0
mpl.rcParams["font.sans-serif"] = "Arial"
mpl.rcParams["savefig.format"] = "pdf"
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import brewer2mpl

#my own
from utils import getlogger, call_sys

#global settings
#logger
date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
logger = getlogger(fn=os.getcwd() + "/" + date.strip() + "_" +
                   os.path.basename(__file__) + ".log")
#data
STAR = "/picb/molsysbio2/caoyaqiang/2.WWX_RNA-seq/0.Tools/STAR-2.5.4b/bin/Linux_x86_64_static/STAR"
STAR_INDEX = "/picb/molsysbio2/caoyaqiang/2.WWX_RNA-seq/1.Reference/3.STAR/150bp"


def concate_files(fnIns, fout):
    logger.info("Merging %s" % (",".join(fnIns)))
    with open(fout, "w") as f:
        for fn in fnIns:
            for line in open(fn):
                f.write(line)
            os.remove(fn)
    logger.info("Merging %s Finished!" % (",".join(fnIns)))


def prepare_sras(metaf=None, SRA_Root=None):
    if metaf == None or SRA_Root == None:
        return
    mat = pd.read_table(metaf, index_col=0, sep="\t")
    data = {}
    for sra in mat.index:
        sample = mat.loc[sra, "LibraryName"]
        sample = sample.split(":")[1]
        sample = sample.replace("CshlLong_RnaSeq_", "")
        sample = sample.replace("CSHL_RnaSeq_", "")
        sample = sample.strip()
        if sample not in data:
            data[sample] = set()
        sra = os.path.join(SRA_Root, sra + ".sra")
        if not os.path.isfile(sra):
            logger.error("For {sample}, {sra} not exists!".format(
                sample=sample, sra=sra))
        data[sample].add(sra)
    for key in data.keys():
        data[key] = list(data[key])
    return data


def dump_sra(sample, sras):
    """
    sample: sras prefix names
    sras: list of sra files
    Only fit for pair-ended sra files. 
    """
    fq1 = sample + "_1.fastq"
    fq2 = sample + "_2.fastq"
    fqs = [fq1, fq2]
    if os.path.exists(fq1) and os.path.exists(fq2):
        logger.info("%s exists!" % ",".join(fqs))
        return fqs
    #dump fastqs
    for i, s in enumerate(sras):
        a = sample + "_" + chr(97 + i)
        cmd = "fastq-dump --split-3 {sra} -A {accession}".format(
            sra=s, accession=a)
        call_sys([cmd])
    fqs = glob.glob("%s*.fastq" % sample)
    groups = {"_1": [], "_2": []}
    for fq in fqs:
        if fq.endswith("_1.fastq"):
            groups["_1"] += [fq]
        if fq.endswith("_2.fastq"):
            groups["_2"] += [fq]
    fqs = []
    for key, value in groups.items():
        fq = sample + key + ".fastq"
        #very important for STAR! Keep the consistant order of paired reads!
        value.sort()
        concate_files(value, fq)
        fqs += [fq]
    return fqs


def backup_fastq(fqs, Fastq_Root="2.Fastq/"):
    nfqs = []
    for fq in fqs:
        nf = os.path.join(Fastq_Root, fq + ".bz2")
        nfqs += [nf]
        if os.path.exists(nf):
            continue
        nnf = fq + ".bz2"
        cmd1 = "bzip2 -z -1 {fq}".format(fq=fq)
        cmd2 = "mv {nnf} {Fastq_Root}".format(nnf=nnf, Fastq_Root=Fastq_Root)
        cmds = [cmd1, cmd2]
        call_sys(cmds)
    return nfqs


def prepare_fastq(
        sample,
        sras,
        SRA_Root="1.SRAs/",
        Fastq_Root="2.Fastq/", ):
    fqs = [
        os.path.join(Fastq_Root, sample + "_1.fastq.bz2"),
        os.path.join(Fastq_Root, sample + "_2.fastq.bz2")
    ]
    flag = (os.path.isfile(fqs[0])) and (os.path.isfile(fqs[1]))
    if flag:
        print sample, " FASTQs has been generated."
    else:
        #convert SRA files to fastq files 
        fqs = dump_sra(sample, sras)
        #compress fastqs 
        fqs = backup_fastq(fqs, Fastq_Root=Fastq_Root)
    fqs.sort()
    return sample, fqs


def prepare_fastq2(Fastq_Root="2.Fastq/", ):
    """
    If the fastq files are well prepared, suitable. 
    """
    fastqs = glob.glob(Fastq_Root + "*.fastq.gz")
    #fastqs = glob.glob( Fastq_Root+"*.fastq" )
    data = {}
    for fq in fastqs:
        s = os.path.split(fq)[1]
        if s.endswith("_R1.fastq.gz"):
            sample = s.replace("_R1.fastq.gz", "")
            if sample not in data:
                data[sample] = [0, 0]
            data[sample][0] = fq
            continue
        if s.endswith("_R2.fastq.gz"):
            sample = s.replace("_R2.fastq.gz", "")
            if sample not in data:
                data[sample] = [0, 0]
            data[sample][1] = fq
            continue
    return data


def sort_index_bam(bam):
    pre = os.path.splitext(bam)[0]
    bai = pre + ".bam.bai"
    c1 = "samtools sort -@ 2 {bam} {pre} ".format(bam=bam, pre=pre)
    c2 = "samtools index {bam} {bai}".format(bam=bam, bai=bai)
    cmds = [c1, c2]
    call_sys(cmds)


def STAR_mapping(sample, fqs, mapping_Root="3.Mapping/"):
    #not finished yet,here mainly using unique mapping
    prefix = os.path.join(mapping_Root, sample, sample + "_")
    bam = prefix + "Aligned.out.bam"
    print bam
    if os.path.isfile(bam):
        logger.info("%s exists!" % bam)
        return
    d = os.path.join(mapping_Root, sample)
    if not os.path.exists(d):
        os.mkdir(d)
    #cmd = "{star} --runMode alignReads --runThreadN 6 --genomeDir {genomeDir} --genomeLoad LoadAndKeep --readFilesIn {fq} --readMatesLengthsIn Equal --readFilesCommand 'zcat -1' --clip5pNbases 10 --clip3pNbases 26 --outFileNamePrefix {prefix} --outSAMtype BAM Unsorted --outSAMstrandField intronMotif --outSAMreadID Number --outFilterMultimapNmax 1 --outFilterMismatchNmax 9 --outFilterMismatchNoverLmax 0.1 --outFilterMatchNminOverLread 0.3 --outFilterScoreMinOverLread 0.3 --alignEndsType EndToEnd".format( star=STAR,genomeDir=STAR_INDEX, fq=" ".join(fqs), prefix=prefix )
    #test for untrimmed bam
    cmd = "{star} --runMode alignReads --runThreadN 6 --genomeDir {genomeDir} --genomeLoad LoadAndKeep --readFilesIn {fq} --readMatesLengthsIn Equal --readFilesCommand 'zcat -1' --clip5pNbases 14 --outFileNamePrefix {prefix} --outSAMtype BAM Unsorted --outSAMstrandField intronMotif --outSAMreadID Number --outFilterMultimapNmax 1 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.1 --outFilterMatchNminOverLread 0.3 --outFilterScoreMinOverLread 0.3 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --alignMatesGapMax 200000 --alignIntronMax 200000 --alignEndsType EndToEnd".format(
        star=STAR, genomeDir=STAR_INDEX, fq=" ".join(fqs), prefix=prefix)
    #test for local
    c = "rm -fvr {pre}_STARtmp".format(pre=prefix)
    #sorting and index bam, bzip unmapped fastq files
    pre = os.path.splitext(bam)[0]
    bai = pre + ".bam.bai"
    samsort = "samtools sort {bam} -T {pre} -o {bam}".format(bam=bam, pre=pre)
    samindex = "samtools index {bam} {bai}".format(bam=bam, bai=bai)
    cmds = [cmd, c, samsort, samindex]
    call_sys(cmds)
    """
    cmds = [  ]
    for m in glob.glob( prefix+"*Unmapped.out*"  ):
        cmds += [ "bzip2 -z -1 {f}".format( f=m )  ]
    call_sys( cmds )
    """


def get_chr_size():
    chrs = {}
    for line in open(CHROM_SIZE):
        line = line.split("\n")[0].split("\t")
        chrs[line[0]] = int(line[1])
    return chrs


def validate_bdg(bdg):
    chrs = get_chr_size()
    nbdg = bdg + ".2"
    with open(nbdg, "w") as f:
        for line in open(bdg):
            line = line.split("\n")[0].split("\t")
            if len(line) < 4:
                continue
            if int(line[1]) >= chrs[line[0]]:
                line[1] = str(chrs[line[0]] - 1)
            if int(line[2]) >= chrs[line[0]]:
                line[2] = str(chrs[line[0]] - 1)
            line = "\t".join(line) + "\n"
            f.write(line)
    cmd = "mv %s %s" % (nbdg, bdg)
    call_sys([cmd])


def bdg2bw(bdg):
    print "Validating :", bdg
    validate_bdg(bdg)
    bw = os.path.splitext(bdg)[0] + ".bw"
    cmd1 = "bedGraphToBigWig {bedgraph} {chromsize} {bigwig}".format(
        bedgraph=bdg, chromsize=CHROM_SIZE, bigwig=bw)
    cmd2 = "rm %s" % bdg
    call_sys([cmd1, cmd2])


def STAR_pileup(bam):
    """
    Pileup reads in BAM format to bedGraph, than convert to bigWig.
    """
    prefix = bam.replace("Aligned.out.bam", "")
    if os.path.exists(prefix + "Signal.Unique.str1.out.bw") and os.path.exists(
            prefix + "Signal.Unique.str2.out.bw"):
        return
    cmd = "{star} --runMode inputAlignmentsFromBAM --runThreadN 5 --genomeDir {genomeDir} --genomeLoad LoadAndKeep --outFileNamePrefix {prefix} --outWigType bedGraph --outWigStrand Unstranded --outWigNorm RPM --inputBAMfile {bam}".format(
        star=STAR, genomeDir=STAR_INDEX, prefix=prefix, bam=bam)
    call_sys([cmd])
    for f in glob.glob(prefix + "Signal.UniqueMultiple.*"):
        os.remove(f)
    for f in glob.glob(prefix + "Signal.Unique.*"):
        bdg2bw(f)


def parse_STAR_log(logs=None, fOut="MappingStat.txt",
                   mapping_root="3.Mapping"):
    suf = "_Log.final.out"
    if logs == None:
        s = mapping_root + "/*/*%s" % suf
        logs = glob.glob(s)
        logs.sort()
    data = {}
    columns = [
        "TotalReads", "UniqueMapRatio", "MultipleMapRatio", "ChimeraMapRatio"
    ]
    for log in logs:
        sample = os.path.basename(log)
        sample = sample.replace(suf, "")
        mat = pd.read_table(log, header=None, sep="\t")
        total_reads = int(mat.iloc[4, 1])
        unique_map_ratio = float(mat.iloc[8, 1].replace("%", ""))
        multiple_map_ratio = float(mat.iloc[23, 1].replace("%", ""))
        chimera_map_ratio = float(mat.iloc[32, 1].replace("%", ""))
        data[sample] = pd.Series(
            [
                total_reads, unique_map_ratio, multiple_map_ratio,
                chimera_map_ratio
            ],
            index=columns)
    data = pd.DataFrame(data)
    data = data.T
    if fOut != None:
        data.to_csv(fOut, sep="\t", index_label="Sample")
    return data


def plot_STAR_stat(data, pre="MappingStat"):
    fig, ax = pylab.subplots(figsize=(11, 4))
    colors = brewer2mpl.get_map('Set3', 'qualitative', 8).mpl_colors
    width = 0.5
    samples = list(data.index)
    ind = np.arange(len(data.index))
    #total reads as bar plot, seperate polyA and NonPolyA
    rs = data.loc[:, "TotalReads"].values
    rs = rs / 1000000.0
    ax.bar(ind, rs, width, color=colors[0])
    ax.set_ylabel("Read Pairs/Million")
    ax.set_xticks(ind + width)
    #ax.set_xticklabels( samples,rotation=90,ha="right" )
    ax.set_xticklabels([])
    for i in ind:
        if i % 2 == 0:
            ax.text(i, 60, samples[i], rotation=90, fontsize=6)
        else:
            ax.text(i, 70, samples[i], rotation=90, fontsize=6)
    #map ratio
    ax2 = ax.twinx()
    umr = data.loc[:, "UniqueMapRatio"].values
    mmr = data.loc[:, "MultipleMapRatio"].values
    ax2.plot(
        ind + width / 2,
        umr,
        color="b",
        marker="s",
        ms=1,
        label="Unique Mapping Ratio ")
    ax2.plot(
        ind + width / 2,
        mmr,
        color="r",
        marker="h",
        ms=1,
        label="Multiple Mapping Ratio")
    for t in ax2.get_yticklabels():
        t.set_color("r")
    ax2.set_ylim([0, 100])
    ax2.set_ylabel("Maping Ratio/%")
    leg = ax2.legend(loc="upper right", fancybox=True)
    #leg.get_frame(  ).set_alpha( 0.5 )
    pylab.savefig(pre + ".pdf", )


def STAR_multiplemapping(sample, fqs, mapping_Root="5.MultipleMapping/"):
    #only suitable for multiple mapping,modified the parameters
    d = os.path.join(mapping_Root, sample)
    if not os.path.exists(d):
        os.mkdir(d)
    else:
        print d
        return
    prefix = os.path.join(mapping_Root, sample, sample + "_")
    bam = prefix + "Aligned.out.bam"
    if os.path.isfile(bam):
        logger.info("%s exists!" % bam)
        return
    #cmd = "{star} --runMode alignReads --runThreadN 5 --genomeDir {genomeDir} --genomeLoad LoadAndKeep --readFilesIn {fq} --readMatesLengthsIn Equal --readFilesCommand 'zcat -1' --outFileNamePrefix {prefix} --outSAMtype BAM Unsorted --outSAMstrandField intronMotif --outSAMreadID Number --outFilterMultimapNmax 20 --alignEndsType EndToEnd --twopassMode Basic --outReadsUnmapped None --chimSegmentMin 12 -chimJunctionOverhangMin 12 --alignSJDBoverhangMin 10 --alignMatesGapMax 10000 --alignIntronMax 100000 --chimSegmentReadGapMax 3 --alignSJstitchMismatchNmax 5 -1 5 5".format(star=STAR, genomeDir=STAR_INDEX, fq=" ".join(fqs), prefix=prefix )
    cmd = "{star} --runMode alignReads --runThreadN 5 --genomeDir {genomeDir} --genomeLoad NoSharedMemory --readFilesIn {fq} --readFilesCommand 'zcat -1' --outFileNamePrefix {prefix} --outSAMtype BAM Unsorted --outSAMstrandField intronMotif --outSAMreadID Number --readMatesLengthsIn NotEqual --outFilterMultimapNmax 20 --alignEndsType Local --twopassMode Basic --outReadsUnmapped None --chimSegmentMin 12 -chimJunctionOverhangMin 12 --alignSJDBoverhangMin 10 --alignMatesGapMax 10000 --alignIntronMax 100000 --chimSegmentReadGapMax 3 --alignSJstitchMismatchNmax 5 -1 5 5 --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3".format(
        star=STAR, genomeDir=STAR_INDEX, fq=" ".join(fqs), prefix=prefix)
    c = "rm -fvr {pre}_STARtmp".format(pre=prefix)
    #sorting and index bam, bzip unmapped fastq files
    pre = os.path.splitext(bam)[0]
    bai = pre + ".bai"
    samsort = "samtools sort {bam} -T {pre} -o {bam}".format(bam=bam, pre=pre)
    samindex = "samtools index {bam} {bai}".format(bam=bam, bai=bai)
    cmds = [cmd, c, samsort, samindex]
    call_sys(cmds)


def trimBed(bed):
    #trim the 4th column of raw bed file to thin it
    n = bed + ".2"
    with open(n, "w") as f:
        for line in open(bed):
            line = line.split("\n")[0].split("\t")
            try:
                line[3] = line[5]
            except:
                report = "FLAG_C an error line in %s: %s" (bed,
                                                           "\t".join(line))
                continue
            line = "\t".join(line) + "\n"
            f.write(line)
    os.system("mv %s %s" % (n, bed))


def bam2bed(bam):
    fd = os.path.splitext(bam)[0]
    bed = fd + ".bed"
    if os.path.exists(bed):
        return
    c = "bamToBed -i {} > {}".format(bam, bed)
    call_sys([c])
    trimBed(bed)


def pipeline():
    mapping_Root = "./"
    #prepare working directories
    Fastq_Root = "../1.Fastqs/"
    data = prepare_fastq2(Fastq_Root)
    #Parallel( n_jobs=8 )( delayed( STAR_mapping )( sample,fqs,mapping_Root  ) for sample,fqs in data.items(  ) )
    Parallel(n_jobs=6)(delayed(STAR_multiplemapping)(sample, fqs, mapping_Root)
                       for sample, fqs in data.items())
    data = parse_STAR_log(
        logs=None, fOut="MappingStat.txt", mapping_root=mapping_Root)
    data = pd.read_table("MappingStat.txt", index_col=0)
    plot_STAR_stat(data, pre="MappingStat")
    #bams = glob.glob( mapping_Root + "*/*out.bam" )
    #bams.extend( glob.glob( multimapping_Root + "*/*out.bam" ) )
    #Parallel( n_jobs=12 )( delayed( STAR_pileup )( bam ) for bam in bams )
    #Parallel( n_jobs=-1 )( delayed( bam2bed )( bam ) for bam in bams )


def main():
    pipeline()


if __name__ == '__main__':
    start_time = datetime.now()
    main()
    elapsed = datetime.now() - start_time
    print "The process is done"
    print "Time used:", elapsed
