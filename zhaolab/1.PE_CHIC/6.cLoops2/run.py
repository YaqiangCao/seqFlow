import os
from glob import glob

cs=["chr%s"%i for i in range(1,23)]
cs.extend(["chrX","chrY"])
cs = ",".join(cs)

for f in glob("../4.reduBedpe/*.bedpe.gz"):
    n = f.split("/")[-1].split(".bedpe.gz")[0]
    if os.path.isdir(n):
        continue
    cmd = "cLoops2 pre -p 10 -f %s -o %s -c %s"%(f,n,cs)
    print(cmd)
    os.system(cmd)

