import gzip
from glob import glob

def bedpe2hicsum(f):
    print f
    with open(f.split("/")[-1]+".hicsum","w") as fo:
        for line in gzip.open(f):
        #for line in open(f):
            line = line.split( "\n" )[ 0 ].split( "\t" ) 
            nline = [line[6],line[0],line[1],line[8],line[3],line[4],line[9]] 
            fo.write("\t".join(nline)+"\n")

#bedpe2hicsum("../../2.bedpesHiCPro/GM12878_HiC_GSM1551552.bedpe.gz")
bedpe2hicsum("../../2.bedpesHiCPro/K562_HiC_GSM1551619.bedpe.gz")
