#!/usr/bin/env python
#--coding:utf-8--
"""
"""

#sys
import os
import time
from glob import glob
from datetime import datetime


import click

@click.command()
@click.option("-name",
              required=True,
              help="Sample name/id for the data.")
@click.option("-org",
              required=True,
              help="Organism for the data.",
              type=click.Choice(["hg38", "mm10"]))
def main(name, org):
    wd = os.path.dirname(os.path.realpath(__file__))  
    print("%s/run.py -org %s"%(wd,org))
    c1 = "cd %s/2.sepSingleDNARNA/"%wd
    c2 = "python coPre.py"
    c3 = "python getStat.py"
    c4 = "cd %s/3.DNA/"%wd
    c5 = "python coDNAPro.py -name %s -org %s"%(name,org)
    c6 = "python getDNAstat.py -name %s"%name
    c7 = "cd %s/4.RNA/"%wd
    c8 = "python coRNAPro.py -name %s -org %s"%(name,org)
    c9 = "python getRNAstat.py -name %s"%name
    c10 = "cd %s/5.statFilter/"%wd
    c11 = "python getSumRaw.py" 
    c12 = "cd %s/6.uniqueDNA/"%wd
    c13 = "python getUniqueDNA.py -name %s"%name
    c14 = "python getBw.py -name %s"%name
    with open("run.sh","w") as fo:
        for c in [c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14]:
            fo.write(c+"\n")
    os.system("bash run.sh")

if __name__ == "__main__":
    startTime = datetime.now()
    main()
    elapsed = datetime.now() - startTime
    print("The process is done.")
    print("Time used:", elapsed)
