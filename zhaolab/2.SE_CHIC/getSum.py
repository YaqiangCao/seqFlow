import pandas as pd

fa = "2.mapping/MappingStat.txt"
fb = "3.beds/bedStat.txt"
fo = "summary.txt"

mata = pd.read_csv(fa,index_col=0,sep="\t")
matb = pd.read_csv(fb,index_col=0,sep="\t")

for c in matb.columns:
    mata[c] = matb[c]

mata.to_csv(fo,sep="\t",index_label="sample")
