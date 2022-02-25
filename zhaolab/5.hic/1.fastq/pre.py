import os
from glob import glob

fs = glob("*R1*.fastq.gz")
for f in fs:
    n = f.split("_R1")[0]
    c1 = "mkdir %s"%n
    c2 = "mv %s %s/"%(f,n)
    c3 = "mv %s_R2.fastq.gz %s/"%(n,n)
    for c in [c1,c2,c3]:
        print(c)
        os.system(c)
