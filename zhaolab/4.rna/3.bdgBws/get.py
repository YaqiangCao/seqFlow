import os
from glob import glob

fs = glob("../2.mapping/*/*_Signal.UniqueMultiple.str1.out.bg")
for f in fs:
    n = f.split("/")[-2]
    if os.path.isfile(n+".bw"):
        continue
    cmd = "ln -s %s %s.bdg"%(f,n)
    print(cmd)
    os.system(cmd)
os.system("python bdg2bw.py -pattern './*.bdg' -org mm10")
os.system("rm *.bdg")
