import os
from glob import glob
fs = glob("../5.tsv/*.tsv.gz")
for f in fs:
    n = f.split("/")[-1].split(".tsv.gz")[0]
    cmd = f"getNucCWOS.py -f {f} -o {n} -cof ~/caoy7/Projects/0.Reference/2.mm10/1.fa/mm10.chrom.sizes -bw &"
    print(cmd)
    os.system(cmd)
