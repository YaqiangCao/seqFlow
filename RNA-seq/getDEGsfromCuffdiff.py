import pandas as pd


def getDegs(f,pre,pcut=1e-2,fccut=1):
    ds = {}
    for i, line in enumerate(open(f)):
        if i == 0:
            continue
        line = line.split( "\n" )[ 0 ].split( "\t" ) 
        if float(line[11]) <= pcut and abs(float(line[9])) >= fccut and float(line[7]) + float(line[8]) >= 1:
            t = "%s|%s|%s"%(line[0],line[2],line[3])
            ds[t] = {line[4]:float(line[7]),line[5]:float(line[8]),"log2(fold_change)":0.0-float(line[9]),"p-value":float(line[11]),"q-value":float(line[12])}
    ds = pd.DataFrame(ds).T
    s = ds["log2(fold_change)"]
    s = s.sort(inplace=0,ascending=1)
    ds = ds.loc[s.index,]
    ds.to_csv("%s_degs.txt"%pre,sep="\t")
    
    koup = [t.split("|")[1] for t in s[s>0].index]
    with open("%s_up.list"%pre,"w") as f:
        f.write("\n".join(koup))
    kodown = [t.split("|")[1] for t in s[s<0].index]
    with open("%s_down.list"%pre,"w") as f:
        f.write("\n".join(kodown))


getDegs("B1_vs_parent/gene_exp.diff","B1")
getDegs("B2_vs_parent/gene_exp.diff","B2")
getDegs("L1_vs_parent/gene_exp.diff","L1")
