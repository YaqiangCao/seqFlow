import os,commands,shutil
from datetime import datetime
import pandas as pd
from glob import glob

def runsys(cmd):
    print cmd
    s = datetime.now()
    _code,output = commands.getstatusoutput(cmd)
    #os.system(cmd)
    e = datetime.now()
    t = e-s
    return t


def runHOMER():
    for f in glob("../1.HiCSummaryFormat/*.hicsum"):
        print f
        fn = f.split("/")[-1].split(".")[0]
        if os.path.exists(fn):
            continue
        cmd1 = "makeTagDirectory {fout} -format HiCsummary {fin}  ".format(fout=fn,fin=f)
        cmd2 =  "findHiCInteractionsByChr.pl {fout} -res 2000 -superRes 10000 -cpu 1 > {fout}.txt".format(fout=fn)
        t1 = runsys(cmd1)
        t2 = runsys(cmd2)
        t = t1 +  t2
        print t


runHOMER()
