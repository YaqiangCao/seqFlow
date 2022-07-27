import gzip
from glob import glob
import pandas as pd
from joblib import Parallel, delayed
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator


def getStat(k):
    ds = {}
    df1 = k + "_DNA_R1.fastq.gz"
    rf1 = k + "_RNA_R1.fastq.gz"
    with gzip.open(df1, "rt") as f:
        for r in FastqGeneralIterator(f):
            rid = r[0]
            cid = "_".join(rid.split("_")[:2])
            if cid not in ds:
                ds[cid] = {
                    "totalDNA": 0,
                    "totalRNA": 0,
                }
            ds[cid]["totalDNA"] += 1
    with gzip.open(rf1, "rt") as f:
        for r in FastqGeneralIterator(f):
            rid = r[0]
            cid = "_".join(rid.split("_")[:2])
            if cid not in ds:
                ds[cid] = {
                    "totalDNA": 0,
                    "totalRNA": 0,
                }
            ds[cid]["totalRNA"] += 1
    return ds


def main():
    fs = glob("*_DNA_R1.fastq.gz")
    ks = [f.split("_DNA_R1.fastq.gz")[0] for f in fs]
    ds = Parallel(n_jobs=min(20, len(ks)))(delayed(getStat)(k) for k in ks)
    data = {}
    for d in ds:
        for k, v in d.items():
            data[k] = v
    data = pd.DataFrame(data).T
    data.to_csv("cellStat.txt", sep="\t")


main()
