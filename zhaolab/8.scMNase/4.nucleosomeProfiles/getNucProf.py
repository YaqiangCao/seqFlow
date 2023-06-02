#!/usr/bin/env python
#--coding:utf-8--
"""
getNucProf.py 
Performs two jobs for pooled single-cells or bulk cells, details listed as following :
1.caculate the nucleosome center-weighted reads occupancy score
2.get the profile of subnucleosome (accessibility signal) 

Get the nucleosome center-weighted reads occupancy score according to:
Brogaard, Kristin, et al. "A map of nucleosome positions in yeast at base-pair resolution." Nature 486.7404 (2012): 496-501.

2023-01-06: basically finished.
2023-01-06: modified the score to normalized to 1, to be better normalized with total reads number; as even read sum up score and odd read sum up score is different; also tunning the fragment size selected to build up the score, lager size, more read coverage, fusszy signals.
2023-04-06: add the subnucleosome 
2023-05-02: add parallel writting to speed up
"""

__author__ = "CAO Yaqiang"
__date__ = "01/06/2023"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#sys
import os
import random
import gzip
import subprocess
from datetime import datetime

#3rd library
import click
import numpy as np
from tqdm import tqdm
from joblib import Parallel, delayed


def isTool(name):
    """
    Check if a tool is on PATH and marked as executable.

    @param name: str, 

    @return: bool, true or false
    """
    from distutils.spawn import find_executable
    if find_executable(name) is not None:
        return True
    else:
        return False


def runCmds(cs):
    """
    Run system commands
    """
    for c in cs:
        print("\t".join([str(datetime.now()), c]))
        output = subprocess.getoutput(c)


def loadChromSizes(f):
    """
    Load chromosome size from UCSC size file. 
    
    @param f: str, file from UCSC such as hg38.chrom.sizes, can be obtained through command fetchChromSizes 
    @return cs: dict, key is chromosome and value is max size
    """
    cs = {}
    for line in open(f):
        line = line.split("\n")[0].split("\t")
        if "alt" in line[0] or "chrUn" in line[0] or "random" in line[
                0] or "chrM" in line[0]:
            continue
        cs[line[0]] = int(line[1])
    return cs


def calcGaussianW(odd=True, hns=73):
    """
    Caculate the Gaussian weight based on the distance to the center and the nucleosome fragment type.

    @param odd: bool, whether the read length is odd or even
    @param hns: int, half nucleosome length to define a nucleosome

    @return: np.array of float, gaussian weights
    """
    ws = []
    if odd:
        for i in range(-hns, hns):
            d = abs(i)
            ws.append(np.exp(-0.5 * d / 20 * d / 20))
    else:
        for i in range(-hns - 1, hns):
            d = abs(i)
            ws.append(0.5 * np.exp(-0.5 * d / 20 * d / 20))
    ws = np.array(ws)
    ws = ws / np.sum(ws)
    return ws


def writeBdg(chrom, cov, fout):
    """
    Write numpy array as bedGraph file.
    
    @param chrom: str, chromosome 
    @param cov: numpy.array of signals
    @param fout: str, bdg file
    """
    with open(fout, "a") as fo:
        i = 0
        while i < len(cov) - 1:
            if cov[i] == 0:  #find the non 0 start
                i += 1
                continue
            #find the same value stop
            j = i + 1
            while j < len(cov) - 1:
                if cov[j] != cov[i]:
                    break
                j += 1
            v = cov[i]
            line = [chrom, i, j, v]
            fo.write("\t".join(list(map(str, line))) + "\n")
            if j == len(cov) - 1:
                break
            i = j


def writeSignals(covs, fout, cpu=1):
    """
    Write {chrom:numpy array} as bedGraph file parallelly.
    """
    report = "\t".join([str(datetime.now()), f"Writing signals to {fout}"])
    print(report)
    tmpwd = str(random.random())
    os.mkdir(tmpwd)
    Parallel(n_jobs=cpu, backend="multiprocessing")(
        delayed(writeBdg)(chrom, cov, tmpwd + "/" + chrom)
        for chrom, cov in tqdm(covs.items()))
    c1 = f"cat {tmpwd}/* > {fout}"
    c2 = f"rm -fvr {tmpwd}"
    runCmds([c1, c2])


def calcCWOS(fin,
             fout,
             csizes,
             normalize="valid",
             cpu=1,
             dfilter=[132, 162],
             hns=73,
             decimals=4):
    """
    Caculating center-weighted reads occupancy scores.

    @param fin: str, tsv file for reads, each record is a nucleosome
    @param csizes: {str:int}, dict of chromosome size to initial the numpy array
    @param normalize: str, valid means only use the number of reads in the length of dfilter, all means use all reads from the file
    @param cpu: int, cpu number for parallel writings
    @param dfilter : list of int, distance filters for each read; fine tuned, do not change; decrease may have better results
    @param hns: int, half nucleosome size to define a nucleosome
    @param decimals: int, effective number of floats to keep after 0
    
    @return covs: {str:np.array} chromsome-size numpy.array for caculated center-weighted occupancy score
    """
    #weight for odd length read
    ow = calcGaussianW(True, hns)
    #weight for even length read
    ew = calcGaussianW(False, hns)
    #occupancy score for coverage
    covs = {}
    #total and valid nucleosomes
    total, valid = 0, 0
    if fin.endswith(".gz"):
        of = gzip.open(fin, "rt")
    else:
        of = open(fin)
    for i, line in enumerate(of):
        if i % 1000000 == 0:
            report = "\t".join(
                [str(datetime.now()), f"{i} records parsed from {fin}"])
            print(report)
        line = line.split("\n")[0].split("\t")
        chrom = line[0]
        start = int(line[1])
        end = int(line[2])
        dis = end - start
        #filter reads not in the chrom.size file
        if chrom not in csizes:
            continue
        total += 1
        if not (dis >= dfilter[0] and dis <= dfilter[1]):
            continue
        valid += 1
        if chrom not in covs:
            covs[chrom] = np.zeros(csizes[chrom], dtype=float)
        center = int((start + end) / 2)
        p = center - hns
        if dis % 2 == 0:
            #even length, two centers
            covs[chrom][p:p + len(ew)] += ew
        else:
            #odd length
            covs[chrom][p:p + len(ow)] += ow
    report = "\t".join([
        str(datetime.now()),
        f"{total} reads and {valid} valid nucleosomes parsed from {fin}."
    ])
    print(report)
    #do the normalization
    report = "\t".join([
        str(datetime.now()), "Normalizing signals and keeping limited decimals"
    ])
    print(report)
    for chrom, cov in tqdm(covs.items()):
        #do the normalization
        if normalize == "valid":
            cov = cov / valid * 10**6
        else:
            cov = cov / total * 10**6
        #only keep limited floats
        cov = cov.round(decimals=decimals)
        covs[chrom] = cov
    writeSignals(covs, fout, cpu=cpu)


def getSubProf(fin,
               fout,
               csizes,
               normalize="valid",
               cpu=1,
               dfilter=80,
               decimals=4):
    """
    Caculating center-weighted reads occupancy scores.

    @param fin: str, tsv file for reads, each record is a nucleosome or subnucleosome
    @param normalize: str, valid means only use the number of reads in the length of dfilter, all means use all reads from the file
    @param dfilter : int, distance filters for each read; fine tuned, do not change; decrease may have better results
    @param decimals: int, effective number of floats to keep after 0
    
    @return covs: {str:np.array} chromsome-size numpy.array
    """
    #occupancy score for subnucleosome coverage
    covs = {}
    #total and valid nucleosomes
    total, valid = 0, 0
    if fin.endswith(".gz"):
        of = gzip.open(fin, "rt")
    else:
        of = open(fin)
    for i, line in enumerate(of):
        if i % 1000000 == 0:
            report = "\t".join(
                [str(datetime.now()), f"{i} records parsed from {fin}"])
            print(report)
        line = line.split("\n")[0].split("\t")
        chrom = line[0]
        start = int(line[1])
        end = int(line[2])
        dis = end - start
        #filter reads not in the chrom.size file
        if chrom not in csizes:
            continue
        total += 1
        if not (dis <= dfilter):
            continue
        valid += 1
        if chrom not in covs:
            covs[chrom] = np.zeros(csizes[chrom], dtype=float)
        covs[chrom][start:end + 1] += 1
        report = "\t".join([
            str(datetime.now()),
            f"{total} reads and {valid} valid subnucleosomes parsed from {fin}."
        ])
    print(report)
    report = "\t".join([
        str(datetime.now()), "Normalizing signals and keeping limited decimals"
    ])
    print(report)
    for chrom, cov in tqdm(covs.items()):
        #do the normalization
        if normalize == "valid":
            cov = cov / valid * 10**6
        else:
            cov = cov / total * 10**6
        #only keep limited floats
        cov = cov.round(decimals=decimals)
        covs[chrom] = cov
    writeSignals(covs, fout, cpu=cpu)


def bdg2bw(bdg, csf):
    """
    Converting .bdg file to .bw file through bedGraphToBigWig.

    @param f: str, .bdg file 
    @param csf: str, chromosome size file
    """
    sbdg = bdg + ".22."
    bw = ".bw".join(bdg.rsplit(".bdg", 1))
    c1 = f"bedSort {bdg} {sbdg}"
    c2 = f"bedGraphToBigWig {sbdg} {csf} {bw}"
    c3 = f"rm {sbdg}"
    runCmds([c1, c2, c3])


@click.command()
@click.option(
    "-f",
    required=True,
    help="Input file in .tsv or .tsv.gz format.",
    type=str,
)
@click.option(
    "-o",
    required=True,
    help="Output file prefix",
    type=str,
)
@click.option(
    "-csf",
    required=True,
    help=
    "Chromosome size file. Can be obtained through the command of fetchChromSizes.",
    type=str,
)
@click.option(
    "-normalize",
    required=False,
    default="valid",
    help=
    "Normalize to total number of reads or valid number of nucleosomes in the size of 132-162 bp. Default is valid.",
    type=click.Choice(["total", "valid"]),
)
@click.option(
    "-bw",
    required=False,
    default=False,
    help=
    "Whether to also output BigWig file. Default is not. If selected, bedSort and bedGraphToBigWig command should be available in the environment. ",
    is_flag=True,
)
@click.option("-p",
              default=4,
              type=int,
              help="Number of CPUs to finish the job, default is set to 4.")
def run(f, o, csf, p=4, normalize="valid", bw=False):
    """
    Caculating the nucleosome center-weighted occupancy score and subnucleosome accessibility from mapped reads.

    Gaussian weights are described in Brogaard, Kristin, et al. "A map of nucleosome positions in yeast at base-pair resolution." Nature 486.7404 (2012): 496-501. The output signals are normalized to total valid reads in the length range of 132-162 bp (+/-15bp).

    Example: 

    getNucProf.py -f test.tsv.gz -o test -cof hg38.chrom.sizes 
    """
    #start
    start = datetime.now()
    cmd = os.path.basename(__file__)
    report = "\t".join([
        str(datetime.now()),
        f"{cmd} -f {f} -o {o} -csf {csf} -normalize {normalize}"
    ])
    print(report)

    #check if output file exists
    foa = o + "_cwos.bdg"
    fob = o + "_subn.bdg"
    if os.path.isfile(foa):
        report = "\t".join(
            [str(datetime.now()), f"Error! Output file {foa} exists! Return."])
        print(report)
        return
    if os.path.isfile(fob):
        report = "\t".join(
            [str(datetime.now()), f"Error! Output file {fob} exists! Return."])
        print(report)
        return

    #caculating CWOS
    csizes = loadChromSizes(csf)
    report = "\t".join([
        str(datetime.now()),
        "Caculating center-weighted nucleosome occupancy score."
    ])
    print(report)
    calcCWOS(f, foa, csizes, normalize=normalize, cpu=p)
    #get the accessibility signal
    report = "\t".join([
        str(datetime.now()),
        "Getting the subnucleosome accessibility profiles."
    ])
    print(report)
    getSubProf(f, fob, csizes, normalize=normalize, cpu=p)

    #file converting
    if bw:
        flag = True
        for t in ["bedSort", "bedGraphToBigWig"]:
            if not isTool(t):
                report = "\t".join([
                    str(datetime.now()),
                    f"Command {t} it not available! Can not convert signal in bedGraph to bigWig file. "
                ])
                print(report)
                flag = False
        if flag:
            report = "\t".join([
                str(datetime.now()), "Converting output file to bigWig file."
            ])
            print(report)
            bdg2bw(foa, csf)
            bdg2bw(fob, csf)
    #finished
    end = datetime.now()
    usedTime = end - start
    report = "\t".join(
        [str(datetime.now()), f"{cmd} job finished. Used time: {usedTime}"])
    print(report + "\n\n")


if __name__ == "__main__":
    run()
