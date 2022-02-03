import os
from glob import glob

fs = glob("../3.allValid/*.gz")
for f in fs:
    n = f.split("/")[-1].split(".allValidPairs")[0]
    cmd = "hicpro2bedpe.py -f %s -o %s &"%(f,n)
    print(cmd)
    os.system(cmd)
