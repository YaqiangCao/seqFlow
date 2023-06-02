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
@click.option("-cpu", default=20, help="CPU numbers to fetch fastqs.")
def main(cpu):
    f = glob("*.txt")[0]
    fqPaths = getFqPaths()
    cmds = []
    ds = pd.read_csv(f, index_col=0, sep="\t")
    for t in ds.itertuples():
        sample = t[0]
        d = ""
        if sample.startswith("GC4654"):
            d = "WT_Thymus-plate1_naiveCD4_"
            sample = "WT_Thymus-plate1_naiveCD4_"+sample
        elif sample.startswith("GC4655"):
            d = "WT_Thymus-plate2_naiveCD4_"
            sample = "WT_Thymus-plate2_naiveCD4_"+sample
        elif sample.startswith("GC4656"):
            d = "WT_Thymus-plate3_naiveCD4_"
            sample = "WT_Thymus-plate3_naiveCD4_"+sample
        elif sample.startswith("GC4657"):
            d = "WT_Thymus-plate4_naiveCD4_"
            sample = "WT_Thymus-plate4_naiveCD4_"+sample
        elif sample.startswith("GC4658"):
            d = "WT_Thymus-plate5_naiveCD4_"
            sample = "WT_Thymus-plate5_naiveCD4_"+sample
        elif sample.startswith("GC4659"):
            d = "WT_Thymus-plate6_naiveCD4_"
            sample = "WT_Thymus-plate6_naiveCD4_"+sample
        elif sample.startswith("GC4660"):
            d = "WT_Thymus-plate7_naiveCD4_"
            sample = "WT_Thymus-plate7_naiveCD4_"+sample
        if t[0] not in fqPaths:
            continue
        if not os.path.exists(d):
            os.mkdir(d)
        fs = fqPaths[t[0]]
        if len(fs) == 1:
            #cmd = "rsync -aP caoy7@137.187.135.165:%s %s.fastq.gz"%(fs[0],gb)
            cmd = "cp %s %s/%s.fastq.gz" % (fs[0],d, sample)
            cmds.append(cmd)
        elif len(fs) == 2:
            cmd1 = "cp %s %s/%s_R1.fastq.gz" % (fs[0],d, sample)
            cmd2 = "cp %s %s/%s_R2.fastq.gz" % (fs[1],d, sample)
            cmds.append(cmd1)
            cmds.append(cmd2)
        else:
            flag = False
            for f in fs:
                if "_R2_" in f:
                    flag = True
            if flag and len(fs)%2 == 0:
                r1s = [f for f in fs if "_R1_" in f]
                cmd1 = "cat %s > %s/%s_R1.fastq.gz"%(" ".join(r1s),d, sample)
                cmd2 = cmd1.replace("_R1","_R2")
                cmds.append(cmd1)
                cmds.append(cmd2)
            else:
                cmd = "cat %s > %s/%s.fastq.gz"%(fa, d,sample)
                cmds.append(cmd)
    Parallel(n_jobs=cpu)(delayed(callSys)(c) for c in cmds)


if __name__ == "__main__":
    main()
