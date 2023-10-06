import os
from glob import glob

ref="/home/caoy7/caoy7/Projects/0.Reference/2.mm10/3.index/2.bowtie2/mm10"
bk="/home/caoy7/caoy7/Projects/0.Reference/2.mm10/11.blacklist/mm10-blacklist.v2.bed"

for d in glob("../1.fastq/*"):
    if os.path.isdir(d):
        n = d.split("/")[-1]
        cmd = f"preIscSeq.py -d {d} -o {n} -barcode barcodes.txt -blacklist {bk} -ref {ref} -p 40"
        print(cmd)
        os.system(cmd)
