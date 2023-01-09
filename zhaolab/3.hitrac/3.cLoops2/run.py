import os
from glob import glob

cs=["chr%s"%i for i in range(1,23)]
cs.extend(["chrX","chrY"])
cs = ",".join(cs)

for f in glob("../2.tracPre2/*/*unique.bedpe.gz"):
    n = f.split("/")[-1].split("_unique.bedpe")[0]
    if os.path.isdir(n):
        continue
    cmd = "cLoops2 pre -p 10 -f %s -o %s -c %s"%(f,n,cs)
    print(cmd)
    os.system(cmd)

"""
for d in glob("*"):
    if os.path.isdir(d):
        if os.path.isfile(d+".hic"):
            continue
        #cmd = "cLoops2 dump -hic -hic_org hg38 -hic_res 200000,50000,20000,5000,1000 -d %s -o %s"%(d,d)
        cmd = "cLoops2 dump -washU -d %s -o %s &"%(d,d)
        print(cmd)
        os.system(cmd)
"""
