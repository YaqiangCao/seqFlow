from glob import glob
import pandas as pd

ds = glob("hicPro/hic_results/stats/*")
data = {}
for d in ds:
    n = d.split("/")[-1]
    i = 0
    data[n] = {}
    for line in open(d+"/"+n+".mpairstat"):
        if line.startswith("#"):
            continue
        line = line.split("\n")[0].split("\t")
        data[n]["%s_%s"%(i, line[0] )] = line[1]
        i += 1
    for line in open(d+"/"+n+"_allValidPairs.mergestat"):
        if line.startswith("#"):
            continue
        line = line.split("\n")[0].split("\t")
        if line[0] == "cis_shortRange":
            line[0] = "cis_shortRange (<20kb)"
        if line[0] == "cis_longRange":
            line[0] == "cis_longRange (>20kb)"
        data[n]["%s_%s"%(i, line[0] )] = line[1]
        i += 1
data = pd.DataFrame(data).T
print(data)
data.to_csv("HiCPro_summary.txt",sep="\t")
