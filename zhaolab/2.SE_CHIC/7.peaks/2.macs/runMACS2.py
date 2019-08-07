#--coding:utf-8--
"""
runMACS2.py
"""

#sys
import os
from glob import glob
from datetime import datetime

#3rd
from joblib import Parallel, delayed


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


def callPeaks(bed):
    sample = bed.split("/")[-1].split(".")[0]
    doMACS = "macs2 callpeak -f BED --nomodel --shift 0 --extsize 75 -t {bed} -n {sample} -B --SPMR -g mm --keep-dup 1 --outdir {out}".format( bed=bed, sample=sample, out="./" + sample)
    #doMACS = "macs2 callpeak -f BEDPE -t {bed} -n {sample} -B --SPMR -g mm --keep-dup 1 --outdir {out}".format( bed=bed, sample=sample, out="./" + sample)
    #doMACS = "macs2 callpeak -f BAMPE --nomodel --broad --shift 0 --extsize 75 -t {bed} -n {sample} -B --SPMR -g mm --keep-dup 1 --outdir {out}".format( bed=bed, sample=sample, out="./" + sample)
    print(doMACS)
    os.system(doMACS)


def main():
    fs = glob("../3.beds/*.bed.gz")
    Parallel(n_jobs=len(fs))(delayed(callPeaks)(f) for f in fs)


if __name__ == "__main__":
    startTime = datetime.now()
    main()
    elapsed = datetime.now() - startTime
    print "The process is done."
    print "Time used:", elapsed
