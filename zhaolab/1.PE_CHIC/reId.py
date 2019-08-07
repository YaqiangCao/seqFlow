import os
from glob import glob

def reid(fs):
    ss = {}
    for line in open("idmap.txt"):
        line = line.split( "\n" )[ 0 ].split( "\t" )
        if len(line) == 1:
            continue
        line[0] = line[0].strip()
        line[1] = line[1].strip()
        #ss[ line[0] ] =line[0]+"_"+ line[1]
        ss[ line[0] ] = line[1]
    for f in fs:
        d = "/".join(f.split("/")[:-1])
        suff = ".".join(f.split("/")[-1].split(".")[1:])
        n = f.split("/")[-1].split(".")[0].split("_")[0]
        if n not in ss:
            continue
        nf = os.path.join(d,ss[n]+"."+suff)
        cmd = "mv %s %s"%(f,nf)
        print(cmd)
        os.system(cmd)
        


#fs = glob("3.1.reduBeds/*.bed.gz")
#fs.extend(glob("5.bdgBws/*.bw"))
fs = glob("5.bdgBws/*.bdg")
fs.sort()
reid(fs)
