import pandas as pd

fa = "1.fastq/5_20190806_KZ1827_fq.txt"
fb = "2.mapping/MappingStat.txt"
fc = "3.beds/bedStat.txt"
fo = "KZ1840_summary.txt"

mata = pd.read_table(fa,index_col=0,sep="\t")
matb = pd.read_table(fb,index_col=0,sep="\t")
matc = pd.read_table(fc,index_col=0,sep="\t")

for c in matb.columns:
    mata[c] = matb[c]
for c in matc.columns:
    mata[c] = matc[c]

mata.to_csv(fo,sep="\t",index_label="sample")
