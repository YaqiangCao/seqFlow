from glob import glob
import pandas as pd

fs = glob("./*/HiCUP_summary_report*.txt")
data = {}
for f in fs:
    print(f)
    mat = pd.read_csv(f,index_col=0,sep="\t")
    #n = mat.index[0].split("/")[-1].split(".hicup.bam")[0]
    n = f.split("/")[-2]
    s = mat.iloc[0,]
    data[n] = s

data = pd.DataFrame( data ) 
cs = list(data.columns)
cs.sort()
data = data[cs].T
data.to_csv("process_stat.txt",sep="\t")
