import os
from glob import glob

fs = glob("../2.pre/*/*_unique.tsv.gz")
for f in fs:
    n = f.split("/")[-2]
    cmd = f"getNucTsvStat.py -f {f} -o {n} -filter -sep -p 10"
    print(cmd)
    os.system(cmd)
    cmd = f"cat {n}_filterSep/*.tsv.gz > {n}.tsv.gz"
    print(cmd)
    os.system(cmd)
