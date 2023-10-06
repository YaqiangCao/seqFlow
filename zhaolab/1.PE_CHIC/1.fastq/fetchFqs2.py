import os
import click
import pandas as pd
from glob import glob
from joblib import Parallel, delayed


def callSys(cmd):
    print(cmd)
    os.system(cmd)


def getFqPaths():
    """
    """
    root = "/mnt/usbdisk/fastq_backup/"
    fs = glob(root + "*/*.fastq.gz")
    ds = {}
    for f in fs:
        n = f.split("/")[-1]
        if n.startswith("KZ"):
            gb = n.split("_")[1]
        else:
            gb = n.split("_")[0]
        #gb = n.split("_")[0]
        if gb not in ds:
            ds[gb] = []
        ds[gb].append(f)
    for gb, sfs in ds.items():
        sfs.sort()
        ds[gb] = sfs
    return ds


@click.command()
@click.option("-f", required=True, help="Samples related meta information.")
@click.option("-cpu", default=50, help="CPU numbers to fetch fastqs.")
def main(f, cpu):
    fqPaths = getFqPaths()
    cmds = []
    ds = pd.read_csv(f, index_col=0, sep="\t")
    for t in ds.itertuples():
        sample = t[0] + "_" + t[11].replace("-", "_").replace(" ", "_")
        #sample = t[0] + "_ChIP-seq_" + t[5]
        #sample = t[12]
        if t[0] not in fqPaths:
            continue
        fs = fqPaths[t[0]]
        if len(fs) == 1:
            #cmd = "rsync -aP caoy7@137.187.135.165:%s %s.fastq.gz"%(fs[0],gb)
            cmd = "cp %s %s.fastq.gz" % (fs[0], sample)
            cmds.append(cmd)
        else:
            flag = False
            for f in fs:
                if "_R2_" in f:
                    flag = True
            if flag and len(fs)%2 == 0:
                r1s = [f for f in fs if "_R1_" in f]
                cmd1 = "cat %s > %s_R1.fastq.gz"%(" ".join(r1s), sample)
                cmd2 = cmd1.replace("_R1","_R2")
                cmds.append(cmd1)
                cmds.append(cmd2)
            else:
                cmd = "cat %s > %s.fastq.gz"%(fa, sample)
                cmds.append(cmd)
    Parallel(n_jobs=cpu)(delayed(callSys)(c) for c in cmds)


if __name__ == "__main__":
    main()
