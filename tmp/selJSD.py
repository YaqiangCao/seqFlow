import glob
import pandas as pd
import numpy as np
import string



def getJSDmat(grps,jsdcut=0.7):
    mat = {}
    for g in grps:
        f = "../1.JSDs/%s.jsd"%g
        print f
        jsd = pd.Series.from_csv(f,sep="\t")
        mat[g] = jsd
    mat = pd.DataFrame(mat)
    ngs = {}
    for t in mat.itertuples():
        i = np.array(t[1:]).argmax()
        if t[i+1] < jsdcut:
            continue
        if mat.columns[i] not in ngs:
            ngs[mat.columns[i]] = []
        ngs[mat.columns[i]].append(t[0])
    ns = []
    for g in grps:
        ns.extend(ngs[g])
    mat = pd.read_table("../../../11.ELRs_PLRs_Binary_Sets/2.ELRs/ELRs.txt",index_col=0) 
    mat = mat.loc[ns,]
    mat.to_csv("selELRs.txt",sep="\t",index_label="TEs")



grps = ["ESC_&_iPSC","ESDR","Blood_&_T_Cell","HSC_&_B_Cell","Brain","Digestive","Epithelia","Mesench","Muscle","Heart","Smooth_Muscle"]
ngs = getJSDmat(grps)
