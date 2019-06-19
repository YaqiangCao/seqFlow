#!/usr/bin/env python2.7
#--coding:utf-8--
"""
"""

__author__ = "CAO Yaqiang"
__date__ = ""
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#sys
import os, json, gzip, time
from glob import glob
from itertools import izip
from datetime import datetime

#3rd
from joblib import Parallel, delayed
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator

seqIds = {
    "PFV": {
        "Pfv89": "GTCAGTCTGTCTTAATACAAAATTCCATGACA",
        "Pfv90": "GGCACATAATGGCGATACAAAATTCCATGACA",
        "Pfv91": "CAACTGGTAAGCAAATACAAAATTCCATGACA",
        "Pfv92": "ACTCTCGCAGTGGGATACAAAATTCCATGACA",
        "Pfv93": "TTAATTCGCCAGGCATACAAAATTCCATGACA",
        "Pfv94": "ATCAGCGTTAACGGATACAAAATTCCATGACA",
        "Pfv95": "CTCTTAATCACGGTATACAAAATTCCATGACA",
        "Pfv96": "GCGACGTACCTCTAATACAAAATTCCATGACA",
    },
    "TNP": {
        "Tnp01": "TGTTCTCCTGAGCAGATGTGTATAAGAGACAG",
        "Tnp02": "ACCTTCGTCCAGATGTGTATAAGAGACAG",
        "Tnp03": "AGAGGTCTCGACAGATGTGTATAAGAGACAG",
        "Tnp04": "GGATGTCAGATGTGTATAAGAGACAG",
        "Tnp05": "TTCCGTGTTAGATGTGTATAAGAGACAG",
        "Tnp06": "ACAACACTCTCACAGATGTGTATAAGAGACAG",
    },
}


def prepare_fastq(Fastq_Root="../2.reid/"):
    """
    If the fastq files are well prepared, suitable. 
    """
    fastqs = glob(Fastq_Root + "*.fastq.gz")
    data = {}
    for fq in fastqs:
        s = os.path.split(fq)[1]
        s = s.replace("_001.fastq.gz", "")
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
        if not s.endswith("_R1") and not s.endswith("_R2"):
            data[s] = [fq]
    for key in data.keys():
        if 0 in data[key] or len(data[key]) != 2:
            del data[key]
    return data


def sepFqs(fq1, fq2, pre):
    """
    Seperate sequences according to PFV and TNP sequences.
    """
    stats = {
        "both": 0,
        "single": {
            "PFV": 0,
            "TNP": 0
        },
        "none": 0,
        "PFV_TNP": {}
    }
    #prepare output files and stats
    fouts = {}
    for p in seqIds["PFV"].keys():
        for t in seqIds["TNP"].keys():
            key = p + "_" + t
            n = "_".join([pre, p, t])
            fouts[key] = {
                "fo_r1": gzip.open(n + "_R1.fastq.gz", "wb"),
                "fo_r2": gzip.open(n + "_R2.fastq.gz", "wb"),
            }
            stats["PFV_TNP"][key] = 0
    #processing pairing fastqs
    #with open(fq1) as f1, open(fq2) as f2:
    with gzip.open(fq1, "rb") as f1, gzip.open(fq2, "rb") as f2:
        i = 0
        #for r1,r2 in zip(SeqIO.parse(f1,"fastq"),SeqIO.parse(f2,"fastq")):
        for r1, r2 in izip(FastqGeneralIterator(f1), FastqGeneralIterator(f2)):
            #r1 = (title,seq,qual)
            r1, r2 = list(r1), list(r2)
            i += 1
            if i % 100000 == 0:
                print("%s reads processed for %s" % (i, pre))
                print(stats)
            pflag, tflag = False, False
            for p, pseq in seqIds["PFV"].items():
                if pseq in r1[1] and r1[1].find(pseq) == 0:
                    pflag = True
                    break
            for t, tseq in seqIds["TNP"].items():
                if tseq in r2[1] and r2[1].find(tseq) == 0:
                    tflag = True
                    break
            #find both id seqs
            if pflag == True and tflag == True:
                ppos = r1[1].find(pseq) + len(pseq) + 1

                tpos = r2[1].find(tseq) + len(tseq) + 1
                r1[1] = r1[1][ppos:]
                r1[2] = r1[2][ppos:]
                r2[1] = r2[1][tpos:]
                r2[2] = r2[2][tpos:]
                stats["both"] += 1
                key = p + "_" + t
                #fouts[key]["counts"] += 1
                fouts[key]["fo_r1"].write("@%s\n%s\n+\n%s\n" %
                                          (r1[0], r1[1], r1[2]))
                fouts[key]["fo_r2"].write("@%s\n%s\n+\n%s\n" %
                                          (r2[0], r2[1], r2[2]))
                stats["PFV_TNP"][key] += 1
            elif pflag == True and tflag == False:
                stats["single"]["PFV"] += 1
            elif tflag == True and pflag == False:
                stats["single"]["TNP"] += 1
            else:
                stats["none"] += 1
    stats["total"] = i
    print(stats)
    with open(pre + "_stat.json", "w") as fo:
        json.dump(stats, fo)


def main():
    """
    sepFqs("test_R1.fastq.gz","test_R2.fastq.gz","test/test")
    """
    data = prepare_fastq("../1.fastq/")
    #specific settings
    ds = {}
    for key, v in data.items():
        if not os.path.exists(key):
            os.mkdir(key)
        nkey = key + "/" + key
        ds[nkey] = v
    data = ds
    Parallel(n_jobs=len(data))(delayed(sepFqs)(fqs[0], fqs[1], sample)
                               for sample, fqs in data.items())


if __name__ == '__main__':
    start_time = datetime.now()
    date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
    main()
    elapsed = datetime.now() - start_time
    print("The process is done")
    print("Time used:", elapsed)
