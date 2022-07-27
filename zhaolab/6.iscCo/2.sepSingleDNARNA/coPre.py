
import os
import gzip
from glob import glob
import pandas as pd
from tqdm import tqdm
from Bio.Seq import Seq
from joblib import Parallel, delayed
from Bio.SeqIO.QualityIO import FastqGeneralIterator


TSEQ = "CCTGCAGG"
barcodes = [ "GGAACTGT", "TCTACTGT", "TCTGTTGT", "TAACGTGT", "CATGGTGT", "AGAGTTGT", "TAAGGTGT", "ATGGCTGT", "TACACTGT", "CAATGTGT", "GATATTGT", "ATGACTGT", "CTAAGTGT", "CCTCATGT", "GACTGTGT", "TCCTGCTG", "CGAGAATG", "CTTCGTTG", "TAGTTCTG", "ATCTGATG", "CACTCTTG", "ACTGCATG", "ATGGACTG", "TCAAGTTG", "CATGTGTG", "CAGTGATG", "CACAGTTG", "GGTCAATG", "CTAGTGTG", "TCATCGTG", "AAGAGGTT", "CATGTATT", "GAGACCTT", "TAGCAGTT", "GTACCGAT", "TTCACGAT", "CAATTCGT", "CATAACTT", "TACAGTGT", "CAAGATAT", "TCCGGTAT", "AAGAATCT", "ATTGGCGT", "ACGTACGT", "GCGTAATT", "AATCGGAT", "TGCAGACA", "ATAGATAC", "CATGTAGA", "AGTTGACC", "ATTAAGCG", "GATGGCTT", "GTCTCCTA", "CGTAATTA", "TCGTCGAT", "ACGTACTC", "CGATTACA", "CATATGCT", "TCAGCTTG", "AGCAATCC", "GATAACCA", "ATGACACC", "CTCTGATT", "GACTAAGA", "AGAATCAG", "TAGACGGA", "GTGAACGT", "TGAGCGAA", "GCTTAGTA", "CAAGTCAC", "TACGCGTT", "TAACCAAG", "CTGCAATC", "GTTATATC", "CCTAGTAG", "TGGACATG", "TGATGCGA", "AGGTTGCT", "GGATCATC", "CAGCTCTT", "AACTGCCA", "CGAACTAC", "ATACGACT", "TCGATTAA", "TCGTTAGC", "TGGCGTAT", "ACACACGT", "TTAAGCAT", "GACGTTAA", "AGGCTTGA", "CTGACGTT", "TTCACTAG", "ATATGCTG", "CAAGGTCA", "GAGCGATA", "AATGACAG",]
barcodes = set(barcodes)
linker="AGAACCATGTCGTCAGTGT"

def pre():
    fs = glob("../1.fastq/*.gz")
    fs.sort()
    ds = {}
    for f in fs:
        n = f.split("/")[-1].split(".fastq.gz")[0]
        n = n.replace("_R1","").replace("_R2","")
        if n not in ds:
            ds[n] = []
        ds[n].append(f)
    return ds


def match(s1, s2, mismatch=2):
    """
    Linker match sequence, allow 2 mismatches.
    """
    s = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            s += 1
    if s > mismatch:
        return False
    else:
        return True


def sep(sample, f1, f2):
    print(sample)
    dfq1 = sample + "_DNA_R1.fastq.gz"
    #dfq2 = sample + "_DNA_R2.fastq.gz"
    rfq1 = sample + "_RNA_R1.fastq.gz"
    #rfq2 = sample + "_RNA_R2.fastq.gz"
    total, dna, rna = 0, 0, 0
    withBarcode = 0
    withLinker = 0
    i = 0
    #with gzip.open(dfq1, "wt") as dfq1, gzip.open(dfq2, "wt") as dfq2, gzip.open(rfq1,"wt") as rfq1, gzip.open(rfq2,"wt") as rfq2:
    with gzip.open(dfq1, "wt") as dfq1, gzip.open(rfq1,"wt") as rfq1:
        with gzip.open(f1, "rt") as f1, gzip.open(f2, "rt") as f2:
            for r1, r2 in zip(FastqGeneralIterator(f1), FastqGeneralIterator(f2)):
                total += 1
                if total % 10000 == 0:
                    print("%s reads parsed from %s"%(total, sample))
                umi = r1[1][:6]
                
                #step 1: check barcode
                #barcode in r2 start position
                b = r2[1][:8]
                if b not in barcodes:
                    continue
                withBarcode += 1

                #step 2: check linker sequences
                tmpa = r2[1][8:8 + len(linker)]
                tmpb = r2[1][9:9 + len(linker)]
                tmpc = r2[1][10:10 + len(linker)]
                tmpd = r2[1][11:11 + len(linker)]
                np = -1
                if match(tmpa, linker):
                    np = 8
                if match(tmpb, linker):
                    np = 9
                if match(tmpc, linker):
                    np = 10
                if match(tmpd, linker):
                    np = 11
                #no find of linker
                if np == -1:
                    continue
                withLinker += 1
                """
                #get varialbel C regions
                p = np + len(linker)
                for j in range(np+len(linker), len(r2[1])):
                    if r2[1][j] != "C":
                        p = j
                        break
                """
                    
                #asiign read number id,sample,barcode cell id,and a number, to make sure unique id
                rid = "_".join([sample, b,umi,str(i)])
                i += 1

                #seperate RNA and DNA
                e = r1[1][6:14]
                #RNA
                if e == TSEQ:
                    rna += 1
                    rfq1.write("@%s\n%s\n+\n%s\n" % (rid, r1[1][14:], r1[2][14:]))
                    #rfq2.write("@%s\n%s\n+\n%s\n" % (rid, r2[1][p:], r1[2][p:]))
                #DNA
                else:
                    dna += 1
                    dfq1.write("@%s\n%s\n+\n%s\n" % (rid, r1[1][6:], r1[2][6:]))
                    #dfq2.write("@%s\n%s\n+\n%s\n" % (rid, r2[1][p:], r1[2][p:]))
    return sample, total, withBarcode, withLinker, dna, rna


def main():
    ds = pre()
    nds = {}
    ds = Parallel(n_jobs=min(20,len(ds)))(delayed(sep)(n, fs[0],fs[1]) for n, fs in ds.items())
    for d in ds:
        n, total, withBarcode, withLinker, dna, rna = d[0],d[1],d[2],d[3],d[4],d[5]
        nds[n] = {
             "total": total,
             "withBarcode":withBarcode,
             "withLinker":withLinker,
             "DNA": dna,
             "RNA": rna,
            }
    nds = pd.DataFrame( nds ).T
    nds.to_csv("summary.txt",sep="\t")


main()
