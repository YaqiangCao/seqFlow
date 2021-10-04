import os
from glob import glob

with open("run.sh","w") as fO:
    ds = glob("../1.fastq/*R1*")
    ds.sort()
    for i,d in enumerate(ds):
        sample = d.split("/")[-1].split("_R1")[0]
        if os.path.exists( sample ):
            continue
        os.mkdir(sample)
        r = open("temp.conf").read()
        r = r.replace("HiC_D345Ctrl_rep1",sample)
        with open(sample+".conf","w") as fo:
            fo.write(r)
        line = "/home/caoy7/caoy7/Packages/HiCUP/hicup --config %s.conf &"%sample
        fO.write(line + "\n")
