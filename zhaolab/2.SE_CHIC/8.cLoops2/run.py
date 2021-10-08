import os
from glob import glob


"""
fs = glob("../7.bedpe/*.bedpe.gz")
for f in fs:
    n = f.split("/")[-1].split(".bedpe")[0]
    cmd = "cLoops2 pre -f %s -o %s &"%(f,n)
    print(cmd)
    os.system(cmd)
"""

ds = [ d for d in glob("*") if os.path.isdir(d)]
for d in ds:
    n = d
    cmd = "cLoops2 callPeaks -d %s -o %s -p 10 -eps 100 -minPts 20 &"%( d,n )
    print(cmd)
    os.system(cmd)
