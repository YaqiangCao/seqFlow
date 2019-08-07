#!/usr/bin/env python
#--coding:utf-8--
"""
bdg2washU.py
Converting .bedgraph file to washU genome browser sorted indexed tracks.
2019-06-28:
"""

#sys
import os, time
from glob import glob
from datetime import datetime

#3rd
import click
from joblib import Parallel, delayed

#seqFlow
from utils import isTool, getLogger, callSys

#global setting
#logger
date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
logger = getLogger(fn=os.getcwd() + "/" + date.strip() + "_" +
                   os.path.basename(__file__) + ".log")




def bdg2washU(f):
    """
    Converting .bdg file to washU track.
    """
    n = f.split("/")[-1].replace(".bdg", "")
    if os.path.isfile(n + "_washU.txt.gz") and os.path.isfile(n+"_washU.txt.gz.tbi"):
        return
    #sort the .bdg file
    cmd1 = "bedSort {bdg} {sbdg}".format(bdg=f, sbdg=n + "_washU.txt")
    #bgzip 
    cmd2 = "bgzip %s"%(n+"_washU.txt")
    #tabix
    cmd3 = "tabix -p bed %s"%(n+"_washU.txt.gz")
    callSys([cmd1,cmd2,cmd3], logger)
    


@click.command()
@click.option(
    "-pattern",
    required=True,
    help="Directory and patterns for the .bg files, for example './mouse*.bdg'."
)
@click.option("-cpu",
              default=10,
              help="Number of CPUs to finish the job, default is set to 10.")
def main(pattern, cpu):
    for t in ["bedSort","bgzip", "tabix"]:
        if not isTool(t):
            logger.error("%s not exits!" % t)
            return
    fs = glob(pattern)
    cpu = min(cpu, len(fs))
    Parallel(n_jobs=cpu)(delayed(bdg2washU)(f) for f in fs)


if __name__ == "__main__":
    startTime = datetime.now()
    main()
    elapsed = datetime.now() - startTime
    print("The process is done.")
    print("Time used:", elapsed)
