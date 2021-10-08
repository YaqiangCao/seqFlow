import os
from glob import glob

fs = glob("../4.reduBeds/*.bed.gz")
for f in fs:
    n = f.split("/")[-1].split(".bed")[0]
    cmd = "getBedpeFBed.py -f %s -o %s &"%(f,n)
    print(cmd)
    os.system(cmd)
