import os
from glob import glob

fs = glob("../2.mapping/*/*_Signal.UniqueMultiple.str1.out.bg")
for f in fs:
    n = f.split("/")[-2]
    cmd = "ln -s %s %s.bdg"%(f,n)
    print(cmd)
    os.system(cmd)
os.system("python bdg2bw.py -pattern './*.bdg' -org mm10")
