#!/usr/bin/env python
#--coding:utf-8--
"""
get10xMatrix.py
Get the matrix from 10x .h5 file and do simple normalization.
"""


#sys
from glob import glob

#3rd
import h5py
import numpy as np
import pandas as pd
from tqdm import tqdm
import scipy.sparse as sp_sparse


def get10xH5Matrix(fh5,mingene=1000,mincell=10,mtratio=0.05,mtprefix="mt-"):
    """
    Get the expression matrix from 10x.h5 file, and do filtering to genes and cells, algo remove cells with too high mt-gene ratio.
    """
    #code from 10x to parse .h5 file
    with h5py.File(fh5,"r") as f:
        if u'version' in f.attrs:
            if f.attrs['version'] > 2:
                raise ValueError('Matrix HDF5 file format version (%d) is an newer version that is not supported by this function.' % version)
        else:
            raise ValueError('Matrix HDF5 file format version (%d) is an older version that is not supported by this function.' % version)
        feature_ids = [x.decode('ascii', 'ignore') for x in f['matrix']['features']['id']]
        feature_names = [x.decode('ascii', 'ignore') for x in f['matrix']['features']['name']]        
        barcodes = list(f['matrix']['barcodes'][:])
        matrix = sp_sparse.csc_matrix((f['matrix']['data'], f['matrix']['indices'], f['matrix']['indptr']), shape=f['matrix']['shape'])
        matrix = matrix.todense()
        gs = [ feature_ids[i] + "|" + feature_names[i] for i in range(len(feature_ids)) ]
        mat = pd.DataFrame( matrix,index=gs,columns=barcodes)
        #first filter mt genes
        mtgenes = [ i for i in mat.index if i.split("|")[1].startswith(mtprefix)]
        ss = mat.sum(axis=0)
        sa = np.sum(mat.loc[mtgenes,],axis=0)
        sr = sa/ss
        sr = sr[ sr< mtratio]
        cs = sr.index
        mat = mat[cs]
        #filter cells
        ns = []
        for c in mat.columns:
            s = mat[c]
            s = s[s> 0]
            if len(s) > mingene:
                ns.append( c)
        mat = mat.loc[:, ns]
        #filter genes
        ns = []
        for t in mat.itertuples():
            if np.sum(t[1:]) > mincell:
                ns.append(t[0])
        mat = mat.loc[ns,]
        return mat



def norm10x(mat, zscore=True, minmax=False):
    """
    Normalization according to https://kb.10xgenomics.com/hc/en-us/articles/115004583806-How-are-the-UMI-counts-normalized-before-PCA-and-differential-expression-
    1) normlize count to median per sample
    2) z-score/minmax for each gene
    """
    ss = np.sum(mat, axis=0)  #total sequencing number of each cell
    ss = np.median(ss, axis=0)
    for c in mat.columns:
        sf = np.sum(mat[c]) / ss
        mat[c] = mat[c] / sf
    #log normalization
    mat = np.log2(mat + 1.0)
    ns = []
    #remove all same genes,maybe all 0.
    for t in mat.itertuples():
        if np.std( t[1:]) > 0.0:
            ns.append( t[0] )
    mat = mat.loc[ns,]
    #gene-wise normalization
    if zscore or minmax:
        for t in tqdm(mat.itertuples()):
            s = np.array(t[1:] )
            if zscore:
                #zscore
                mat.loc[t[0],] =  (s - s.mean()) / s.std()
            if minmax:
                #normalized to 0-1
                mat.loc[t[0],] =  (s-s.min())/(s.max()-s.min())
    return mat  #each row a gene, each column a cell

def main():
    for f in glob("../2.count/*/outs/filtered_feature_bc_matrix.h5")
        print(f)
        n = f.split("/")[-3]
        mat = get10xH5Matrix( f )
        mat = norm10x( mat,False,False)
        print(mat.shape)
        mat.to_csv(n+".txt")

main()
