import pylab
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm


def parseTarget(f="DHS_annot_gene.txt",cut=2000):
    r2g,g2r = {},{}    #region to gene or gene to region
    for line in open(f):
        line = line.split("\n")[0].split("\t")
        t = line[1]
        t = t.split("_")
        if len(t) !=5:
            continue
        r = "_".join( t[:3] )
        d = int(t[3])
        g = t[4]
        if abs( d ) > cut:
            continue
        if r not in r2g:
            r2g[r] = []
        if g not in g2r:
            g2r[g] = []
        r2g[r].append( g )
        g2r[g].append( r )
    #print( len(g2r), len(r2g))
    return r2g,g2r



def getExp(f="scrna_cd4th_from_monocle3.txt",cut=5,librarySize=2000):
    mat = {}
    j = 0
    for i,line in enumerate(open(f)):
        if i == 0:
            ss = line.split("\n")[0].split("\t")
            for s in ss:
                mat[s] = {}
        else:
            line = line.split("\n")[0].split("\t")
            g = line[0]
            vs = list( map( float,line[1:])  )
            if np.sum(vs) < cut:
                continue
            j += 1
            for i, v in enumerate( vs ):
                s = ss[i]
                mat[ ss[i] ][ g ] = v
    for s in mat.keys():
        ss = np.array( mat[s].values() )
        if np.sum( ss ) < librarySize:
            del mat[s]
    mat = pd.DataFrame(mat)
    mat.to_csv( "exp.txt",sep="\t")


def sortExp( f="exp.txt",bins=100):
    mat = pd.read_csv( f, sep="\t",index_col=0)
    #print(mat.shape)
    for c in mat.columns:
        mat[c] = mat[c] / mat[c].sum()
    s = mat.std( axis = 1) /  mat.mean( axis=1)
    s = s.sort_values(ascending=False)
    s = s[s>0]
    #print(s[:5])
    gs = {}
    step = len(s) / bins
    for i in range(0,bins-1):
        ts = s.index[ i * step: (i+1)*step ]
        gs[i] =ts
    return gs


def getAcMat(f,fs,sampling=False):
    """
    mat = pd.read_csv(f,index_col=0,sep="\t")
    """
    mat = {}
    for i,line in enumerate(open(f)):
        #if i % 10000 == 0:
        #    print(i)
        if i == 0:
            ss = line.split("\n")[0].split("\t")
            for s in ss:
                mat[s] = {}
        else:
            line = line.split("\n")[0].split("\t")
            g = line[0]
            vs = list( map( int,line[1:])  )
            if np.sum(vs) == 0:
                continue
            for i, v in enumerate( vs ):
                s = ss[i]
                mat[ ss[i] ][ g ] = v
    mat = pd.DataFrame(mat)
    ss = {}
    for line in open(fs):
        line = line.split("\n")[0].split("\t")
        if len(line) !=2:
            continue
        ss[line[1]] = int( line[0] )
    for c in mat.columns:
        mat[c] = mat[c] / ss[c]
    if sampling:
        cs = mat.columns
        np.random.shuffle(cs)
        cs = cs[:46] 
        mat = mat[cs]
    return mat



def getCVcor(g2r,gs,acMat,f="exp.txt"):
    expMat = pd.read_csv(f,index_col=0,sep="\t")
    for c in expMat.columns:
        expMat[c] = expMat[c] / expMat[c].sum()
    
    expcv = expMat.std( axis = 1 ) / expMat.mean( axis=1)
    accv = acMat.std( axis = 1 ) / acMat.mean( axis= 1)
    
    #accv = accv.fillna(0.0)
    #accv = accv[accv>0]
   
    expcv.to_csv("test_expcv.txt",sep="\t")
    accv.to_csv("test_accv.txt",sep="\t")

    x,y = [],[] #x, accessibility; y, expression

    for i, g in gs.items():
        nexpcv = expcv[ g ]
        rs = []
        for t in g:
            if t in g2r :
                for r in g2r[t]:
                    if r in accv.index:
                        rs.extend( g2r[t] )
        naccv = accv[ rs ]
        #print(nexpcv.mean(), naccv.mean())
        x.append(  naccv.mean() )
        y.append(  nexpcv.mean() ) 
    print( np.corrcoef( x,y)[0][1] )
    return x,y



#getExp() #prepare expression matrix
r2g,g2r = parseTarget()


gs = sortExp()
mat = getAcMat( "ataccd4_rc_mat.txt", "wc_uniq_CD4.txt")
getCVcor( g2r, gs, mat )

mat = getAcMat( "clusTcell_rc_mat.txt", "wc_uniq_Tcell.txt")
print(mat.shape)
getCVcor( g2r, gs, mat )
