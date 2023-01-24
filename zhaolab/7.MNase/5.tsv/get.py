import os
import gzip
from glob import glob

fs = glob("../4.reduBedpe/*.bedpe.gz")
for f in fs:
    nf = f.split("/")[-1].split(".bedpe")[0]
    with gzip.open(nf+".tsv.gz","wt") as fo:
        for line in gzip.open(f,"rt"):
            line = line.split("\n")[0].split("\t")
            nline = [line[0],line[1],line[5]]
            fo.write("\t".join(nline)+"\n")
