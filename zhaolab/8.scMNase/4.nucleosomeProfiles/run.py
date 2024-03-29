import os
from glob import glob
from joblib import Parallel, delayed


def runCmds(cs):
    for c in cs:
        print(c)
        os.system(c)


if __name__ == "__main__":
    ds = {
        "DHSs":
        "/home/caoy7/caoy7/Projects/22.NaivePrimedCD4/2.processes/0.referenceData/3.DNase/4.peaks/DNase_Naive_peaks_summit.txt",
        "CTCF":
        "/home/caoy7/caoy7/Projects/22.NaivePrimedCD4/2.processes/0.referenceData/2.CTCF/3.updatedPeaks/WT_Naive_CTCF_peaks_summit.txt",
        "TSS":
        "/home/caoy7/caoy7/Projects/0.Reference/2.mm10/2.annotations/gencode_vM21_pcRNA_tss.txt",
        "immuneCellAllDHSs":
        "/export/1/caoy7_7T/22.NaivePrimedCD4/1.pub/3.ImmuneCells_ATAC/0.peaks_from_paper/ImmGenATAC_mm10_submit.bed"
    }
 
    cmds = []
    fs = glob("../3.tsvStat/*.tsv.gz")
    for f in fs:
        n = f.split("/")[-1].split(".tsv")[0]
        if os.path.isfile(n+"_cwos.bw") and os.path.isfile(n+"_subn.bw"):
            continue
        cs = []
        cmd = f"getNucProf.py -f {f} -o {n} -csf ~/caoy7/Projects/0.Reference/2.mm10/1.fa/mm10.chrom.sizes -bw -p 1"
        cs.append(cmd)
        for k, cf in ds.items():
            if k == "TSS":
                cmd = f"plotNucProf.py -c {cf} -n {n}_cwos.bw -s {n}_subn.bw -o {n}_{k} -title '{n} {k}' -ext 1000 -d"
            else:
                cmd = f"plotNucProf.py -c {cf} -n {n}_cwos.bw -s {n}_subn.bw -o {n}_{k} -title '{n} {k}' -ext 1000 "
            cs.append(cmd)
        cs.append(cmd)
        cmds.append(cs)
    Parallel(n_jobs=6, backend="multiprocessing",
             temp_folder="./tmp/")(delayed(runCmds)(cs) for cs in cmds)
