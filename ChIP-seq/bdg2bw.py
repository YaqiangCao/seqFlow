#--coding:utf-8--
"""
bdg2bw.py
2019-05-13:
"""

#sys
import os
from glob import glob
from datetime import datetime

#3rd
from tqdm import tqdm
from joblib import Parallel, delayed

CHROM = "/mnt/data/tangq/Projects/0.Reference/1.hg38/1.fa/hg38.chrom.sizes"


def callSys(cmds):
    """
    Call systematic commands without return.
    """
    for c in cmds:
        print(c)
        try:
            os.system(c)
        except:
            print("ERROR!")


def getChrSize():
    chrs = {}
    for line in open(CHROM):
        line = line.split("\n")[0].split("\t")
        chrs[line[0]] = int(line[1])
    return chrs


def validateBdg(bdg):
    chrs = getChrSize()
    nbdg = bdg + ".2"
    with open(nbdg, "w") as f:
        for line in open(bdg):
            line = line.split("\n")[0].split("\t")
            if len(line) < 4:
                continue
            if line[0] not in chrs:
                continue
            if int(line[1]) >= chrs[line[0]] or int(line[2]) > chrs[line[0]]:
                continue
            line = "\t".join(line) + "\n"
            f.write(line)
    cmd = "mv %s %s" % (nbdg, bdg)
    os.system(cmd)


def bdg2bw(f):
    n = f.split("/")[-2]
    cmd1 = "bedSort {bdg} {sbdg}".format(bdg=f, sbdg=n + ".bdg")
    callSys([cmd1])
    #validation bdg
    validateBdg(n + ".bdg")
    cmd2 = "bedGraphToBigWig {bdg} {chrom} {bw}".format(bdg=n + ".bdg",
                                                        chrom=CHROM,
                                                        bw=n + ".bw")
    callSys([cmd2, "rm %s.bdg" % n])


def main():
    fs = glob("../2.peaks/*/*.bdg")
    Parallel(n_jobs=len(fs))(delayed(bdg2bw)(f) for f in fs)


if __name__ == "__main__":
    startTime = datetime.now()
    main()
    elapsed = datetime.now() - startTime
    print("The process is done.")
    print("Time used:", elapsed)
