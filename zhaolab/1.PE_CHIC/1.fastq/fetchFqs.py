import os
import click
import pandas as pd
from glob import glob
from joblib import Parallel, delayed


def callSys(cmd):
    print(cmd)
    os.system(cmd)

@click.command()
@click.option("-f",required=True,help="Samples related meta information, contains the fastq file location.")
@click.option("-cpu",default=10,help="CPU numbers to fetch fastqs.")
def main(f,cpu):
    cmds = []
    ds = pd.read_csv(f,index_col=0,sep="\t")
    for t in ds.itertuples():
        gb = t[0]+"_"+t[9].replace("-","_")
        fs = t[-1].split(",")
        if len(fs) == 1:
            cmd = "rsync -aP caoy7@137.187.135.165:%s %s.fastq.gz"%(fs[0],gb)
            cmds.append(cmd)
        elif len(fs) == 2:
            cmd1 = "rsync -aP caoy7@137.187.135.165:%s %s_R1.fastq.gz"%(fs[0],gb)
            cmd2 = "rsync -aP caoy7@137.187.135.165:%s %s_R2.fastq.gz"%(fs[1],gb)
            cmds.append(cmd1)
            cmds.append(cmd2)
        else:
            print("ERROR. Files for %s is %s"%(gb,fs))
    Parallel(n_jobs=cpu)(delayed(callSys)(c) for c in cmds)
    

if __name__ == "__main__":
    main()
