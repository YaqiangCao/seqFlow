import click
import pandas as pd
from glob import glob

def getLibFqs(fs):
    """
    Prepare fastq files into GB1234: [GB1234_R1.fastq.gz,GB1234_R2.fastq.gz]
    """
    fs.sort()
    ds = {}
    for f in fs:
        gb = f.split("/")[-1].split("_")[0]
        if gb not in ds:
            ds[gb] = []
        ds[gb].append( f )
    for gb,sfs in ds.items():
        sfs.sort()
        ds[gb] = sfs
    return ds

def findFqs(meta):
    """
    Find all related fastqs according to specific GB/GC library ID.
    """
    root = "/mnt/usbdisk/fastq_backup/"
    fs = glob(root+"*/*.fastq.gz")
    ds = getLibFqs(fs)
    s = {}
    for g in meta.index:
        fqs = ds[g]
        if len(fqs) == 0:
            print("ERROR!No fastqs were found for %s"%g)
        s[g] = ",".join(fqs)
    meta["fastqs"] = pd.Series(s)
    return meta
    

@click.command()
@click.option("-f",required=True,help="Sample related txt file. Columns should be exact as following:\nLibraryID\tSampleID\tLibraryType\tReferenceGenome\tSampleDiscription")
def main(f):
    meta = pd.read_csv(f,index_col=0,sep="\t")
    meta = findFqs(meta)
    meta.to_csv( f.replace(".txt","_fq.txt"),sep="\t")

if __name__ == "__main__":
    main()
