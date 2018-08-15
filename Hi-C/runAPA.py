import os,commands,shutil
from glob import glob
from datetime import datetime
import pandas as pd
import matplotlib as mpl
mpl.use("pdf")
mpl.rcParams["pdf.fonttype"] = 42
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import pylab
#import brewer2mpl
import seaborn as sns
sns.set_style("white")


def runsys(cmd):
    print cmd
    s = datetime.now()
    #_code,output = commands.getstatusoutput(cmd)
    os.system(cmd)
    e = datetime.now()
    t = e-s
    return t

def runAPA():
    fs = glob("../../1.sets/*juicebox.txt")
    fs.sort()
    for f in fs:
        n = f.split("/")[-1].split(".")[0].split("_")
        tool = n[-1]
        if os.path.exists(tool):
            continue
        hic = "../../0.HiC/GM12878.hic"
        cmd = "java -jar /picb/molsysbio/usr/caoyaqiang/4.ENV/SG/sage-7.3-cyq/upstream/bio/juice/juicer_tools_0.7.5.jar apa -n 0 -w 5 -r 5000 -u {hic} {loop} {output}".format(hic=hic,loop=f,output=tool)
        runsys(cmd)


def statAPAP2L():
    ds = {}
    i = 0
    fs = glob("*/5000/gw/measures.txt")
    #fs.sort()
    for f in fs:
        #n = f.split("/")[0].split("_")
        n = f.split("/")[0]
        print n
        s = pd.Series.from_csv(f,sep="\t")["P2LL"]
        #ds[i] = {"tech":n[0],"tool":n[1],"P2LL":s}
        ds[i] = {"tool":n,"P2LL":s}
        i += 1
    ds = pd.DataFrame(ds).T
    fig,ax = pylab.subplots(figsize=(4,2.75))
    #sns.barplot(x="tool",y="P2LL",hue="tech",data=ds,saturation=0.5,palette="Set3")
    sns.barplot(x="tool",y="P2LL",data=ds,saturation=0.5,palette="Set3")
    pylab.savefig("all_P2LL")


def plotcenterNormedAPA():
    fs = glob("*/5000/gw/centerNormedAPA.txt")
    fs.sort()
    for f in fs:
        n = f.split("/")[0]
        mat = []
        for line in open(f):
            line = line.split( "[" )[ 1 ].split( "]" )[0].split(",")
            mat.append(map(float,line))
        es = open("%s/5000/gw/enhancement.txt"%n).read().split("\n")[0].split()
        es = np.array(map(float,es))
        es[np.isnan(es)] = 0.0
        es[np.isinf(es)] = np.max(es[~np.isinf(es)])
        #es = es[~np.isinf(es)]
        mat = np.array(mat)
        fig,ax = pylab.subplots(figsize=(3,2))
        sns.heatmap(mat,xticklabels=False,yticklabels=False,square=True)
        ax.set_title("%s loops, P2M:%.3f"%(len(es),es.mean()),fontsize=8)
        pylab.savefig("%s.pdf"%n)
        print n, es.mean()

runAPA()
plotcenterNormedAPA()
#statAPAP2L()
