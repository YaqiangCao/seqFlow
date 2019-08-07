import pandas as pd

fa = "1.fastq/4_20190719_KZ1854_fq.txt"
fb = "2.mapping/MappingStat.txt"
fc = "3.bedpe/bedpeStat.txt"
fo = "KZ1854_summary.txt"

mata = pd.read_csv(fa,index_col=0,sep="\t")
matb = pd.read_csv(fb,index_col=0,sep="\t")
matc = pd.read_csv(fc,index_col=0,sep="\t")

for c in matb.columns:
    mata[c] = matb[c]
for c in matc.columns:
    mata[c] = matc[c]

mata.to_csv(fo,sep="\t",index_label="sample")
