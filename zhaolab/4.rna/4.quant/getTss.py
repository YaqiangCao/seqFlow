import pandas as pd
gtf = "/home/caoy7/caoy7/Projects/0.Reference/2.mm10/2.annotations/gencode_vM21_pcRNA_tss.bed" #this file already extend 2500bp up and down
mat = pd.read_csv("KI_Rx1KI_Exp_filter.txt",index_col=0,sep="\t")
gs = set( [ "|".join(t.split("|")[:-1]) for t in  mat.index] )

ext = 1000
with open("gene_tss.bed","w") as fo:
    for line in open(gtf):
        line = line.split("\n")[0].split("\t")
        g = line[3]
        g = g.split("|")
        g[0] = g[0].split(".")[0]
        g = "|".join(g)
        if g in gs:
            line[3] = g
            m = int((int(line[1]) + int( line[2] )) / 2)
            line[1] = str( m - ext )
            line[2] = str( m + ext )
            fo.write( "\t".join(line) + "\n")
