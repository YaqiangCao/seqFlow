import os
from glob import glob
import pandas as pd

os.system("md5sum *.bw *.gz > md5.txt")
md5 = {}
for line in open("md5.txt"):
    line = line.split("\n")[0].split()
    md5[ line[1] ] = line[0]

samplesTable = {}
processDataTable = {}
rawFileTable = {}
pairedTable = {}

fs = glob("*.bw")
for f in fs:
    n = f.split(".bw")[0]
    rf = n+".fastq.gz"
    source = "Naive CD4+ T cells"
    samplesTable[n] = {
            "1_title": n,
            "2_source_name": source,
            "3_organism": "Mus musculus",
            "4_characteristics: strain":"mixed C57BL/6 and 129 ",
            "5_characteristics: cell type":source,
            "6_characteristics: tissue":source,
            "7_characteristics: genotype":"",
            "8_molecule":"DNA",
            "9_description":"ChIP-seq",
            "10_processed_data_file":f,
            #"11_raw_file":r1,
            #"12_raw_file":r2,
            }
    processDataTable[ n ] = {
                "1_file_name": f,
                "2_file_type": "BIGWIG",
                "3_file_checksum": md5[ f ],
            }
    #single end
    if os.path.isfile(rf):
        samplesTable[n]["11_raw_file"] = rf
        rawFileTable[ rf ] = {
                "1_file_name": rf, 
                "2_file_type": "FASTQ",
                "3_file_checksum": md5[ rf ],
            }
    #paired-end
    else:
        r1 = n + "_R1.fastq.gz"
        r2 = n + "_R2.fastq.gz"
        if not os.path.isfile( r1 ) or not os.path.isfile( r2 ):
            print(n,"no fastq files")
            continue
        samplesTable[n]["11_raw_file"]=r1
        samplesTable[n]["12_raw_file"]=r2
        rawFileTable[ r1 ] = {
                "1_file_name": r1, 
                "2_file_type": "FASTQ",
                "3_file_checksum": md5[ r1 ],
            }
        rawFileTable[ r2 ] = {
                "1_file_name": r2, 
                "2_file_type": "FASTQ",
                "3_file_checksum": md5[ r2 ],
            }
        pairedTable[ n ] = {
                "1_file_name_1": r1,
                "2_file_name_2": r2,
            }


samplesTable = pd.DataFrame( samplesTable ).T
samplesTable.to_csv("1_samplesTable.txt",sep="\t")
print(samplesTable)
processDataTable = pd.DataFrame( processDataTable ).T
processDataTable.to_csv("2_processDataTable.txt",sep="\t")
print(processDataTable)
rawFileTable = pd.DataFrame( rawFileTable ).T
rawFileTable.to_csv("3_rawFileTable.txt",sep="\t")
print(rawFileTable)
pairedTable = pd.DataFrame( pairedTable ).T
pairedTable.to_csv("4_pairedTable.txt",sep="\t")
print(pairedTable)
