import os
import gzip
from glob import glob
from collections import Counter
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import matplotlib
matplotlib.use('pdf')
import pylab
import pandas as pd
import seaborn as sns
from joblib import Parallel, delayed


def preFqs(fastqRoot):
    """
    If the fastq files are well prepared, suitable. 
    """
    fastqs = glob(fastqRoot + "/*.fastq.gz")
    data = {}
    for fq in fastqs:
        s = os.path.split(fq)[1]
        s = s.replace(".fastq.gz", "")
        if s.endswith("_R1"):
            sample = s.replace("_R1", "")
            if sample not in data:
                data[sample] = [0, 0]
            data[sample][0] = fq
        if s.endswith("_R2"):
            sample = s.replace("_R2", "")
            if sample not in data:
                data[sample] = [0, 0]
            data[sample][1] = fq
    for key, fqs in data.items():
        if len(fqs) != 2:
            logger.error(
                "for %s there is not paired fastq files, only %s found" %
                (key, ",".join(fqs)))
            del data[key]
    return data


def findLinker(seq, linker):
    """
    Match the linker in the read sequence.
    """
    pos = -1
    for i in range(len(seq) - 9):
        seed = seq[i:i + 9]
        if linker.startswith(seed):
            pos = i
            break
    return pos


def checkLinkerStart(fq1,
                     fq2,
                     pre,
                     linker="CTGTCTCTTATACACATCT",
                     starts=["AATT", "CATG"]):
    """
    Cut linkers and filter too short reads
    """
    p1s = []
    p2s = []
    tot = 0
    l1s = 0
    l2s = 0
    l12s = 0
    r1s = 0  #right starts
    r2s = 0
    r12s = 0
    with gzip.open(fq1, "rt") as f1, gzip.open(fq2, "rt") as f2:
        for r1, r2 in zip(FastqGeneralIterator(f1), FastqGeneralIterator(f2)):
            r1, r2 = list(r1), list(r2)
            tot += 1
            #check starts
            flag1, flag2 = False, False
            if r1[1].startswith(starts[0]) or r1[1].startswith(starts[1]):
                flag1 = True
                r1s += 1
            if r2[1].startswith(starts[0]) or r2[1].startswith(starts[1]):
                flag2 = True
                r2s += 1
            if flag1 and flag2:
                r12s += 1
            #check the linker
            r1pos = findLinker(r1[1], linker)
            r2pos = findLinker(r2[1], linker)
            #trim reads
            if r1pos != -1:
                p1s.append(r1pos)
                l1s += 1
            if r2pos != -1:
                p2s.append(r1pos)
                l2s += 1
            if r1pos != -1 and r2pos != -1:
                l12s += 1
    fig, ax = pylab.subplots()
    p1s = pd.Series(Counter(p1s))
    p2s = pd.Series(Counter(p2s))
    #sns.kdeplot(p1s, label="R1 linker position",ax=ax)
    #sns.kdeplot(p2s, label="R2 linker position",ax=ax)
    #ax.set_xlim([0,50])
    ax.scatter(p1s.index, p1s.values, color="red", label="R1 linker")
    ax.scatter(p2s.index, p2s.values, color="blue", label="R2 linker")
    ax.set_yscale("log")
    ax.set_ylabel("read counts")
    ax.legend()
    ax.set_title(
        "total:%s;R1 linker:%.3f; R2 liner:%.3f; both:%.3f\n R1 right starts: %.3f;R2 right starts:%.3f; both:%.3f"
        % (tot, float(l1s) / tot, float(l2s) / tot, float(l12s) / tot,
           float(r1s) / tot, float(r2s) / tot, float(r12s) / tot))
    pylab.savefig(pre + "_linkers.pdf")
    s = {
        "total PETs": tot,
        "R1 with linker": l1s,
        "R2 with linker": l2s,
        "R1 linker ratio": float(l1s) / tot,
        "R2 linker ratio": float(l2s) / tot,
        "both R1 and R2 linker ratio": float(r12s) / tot,
        "R1 with right start": r1s,
        "R2 with right start": r2s,
        "R1 with right start ratio": float(r1s) / tot,
        "R2 with right start ratio": float(r2s) / tot,
        "both R1 and R2 with right start ratio": float(r12s) / tot,
    }
    return pre, s


ds = preFqs("../1.fastq")
data = Parallel(n_jobs=10)(delayed(checkLinkerStart)(fs[0], fs[1], sample)
                           for sample, fs in ds.items())
ds = {}
for d in data:
    ds[d[0]] = d[1]
ds = pd.DataFrame(ds).T
ns = list(ds.index)
ns.sort()
ds = ds.loc[ns, ]
ds.to_csv("LinkerStart_stat.txt", sep="\t")
