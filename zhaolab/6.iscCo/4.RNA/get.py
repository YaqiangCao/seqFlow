import os
from glob import glob

key = "RNA"
fs = glob("../2.sepSingleDNARNA/*%s_R1.fastq.gz"%key)
fs.sort()
cmd = "cat %s > %s_R1.fastq.gz"%(" ".join(fs),key)
print(cmd)
os.system(cmd)
cmd = cmd.replace("_R1.fastq.gz","_R2.fastq.gz")
print(cmd)
os.system(cmd)
