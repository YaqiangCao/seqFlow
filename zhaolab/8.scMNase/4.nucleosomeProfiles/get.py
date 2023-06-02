import os
from glob import glob

fs = glob("../2.pre/*/*unique.tsv.gz")
ns = [ f.split("/")[-2] for f in fs]
for f in fs:
    n = f.split("/")[-1].split(".tsv.gz")[0]
    cmd = f"python3 getNucProf.py -f {f} -o {n} -csf ~/caoy7/Projects/0.Reference/2.mm10/1.fa/mm10.chrom.sizes.main -bw -p 4 "
    print(cmd)
    os.system(cmd)


ds = {"DHSs":"/home/caoy7/caoy7/Projects/22.NaivePrimedCD4/2.processes/0.referenceData/3.DNase/4.peaks/DNase_Naive_peaks_summit.txt",
"CTCF":"/home/caoy7/caoy7/Projects/22.NaivePrimedCD4/2.processes/0.referenceData/2.CTCF/3.updatedPeaks/WT_Naive_CTCF_peaks_summit.txt"
}

for pre in ns:
    for k, f in ds.items():
        cmd = f"python plotNucProf.py -c {f} -n {pre}_unique_cwos.bw -s {pre}_unique_subn.bw -o {pre}_{k} -title '{pre} {k}' -ext 1000 &"
        os.system(cmd)
