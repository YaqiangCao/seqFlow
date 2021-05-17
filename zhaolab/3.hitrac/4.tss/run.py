import os
from glob import glob

ref = "/home/caoy7/caoy7/Projects/0.Reference/1.hg38/2.annotations/gencode_v30_pcRNA_tss.bed"
#ref = "/home/caoy7/caoy7/Projects/0.Reference/2.mm10/2.annotations/gencode_vM21_pcRNA_tss.bed"
for d in glob("../3.cLoops2/*"):
    if not os.path.isdir(d):
        continue
    n = d.split("/")[-1]
    cmd = "cLoops2 agg -p 20 -d %s -peaks %s -o %s_tss -peak_ext 2500 -peak_bins 200 -peak_norm -skipZeros"%(d,ref,n)
    print(cmd)
    os.system(cmd)
