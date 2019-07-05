import os
from glob import glob
import pandas as pd
import numpy as np


def txt2cdt(f):
    """
    Converting matrix txt file to CDT file that can be loaded in Treeview to generate heatmaps.
    """
    nf = f.replace(".txt", ".cdt")
    if os.path.isfile(nf):
        return
    mat = pd.read_table(f, index_col=0)
    genes, mat, columns = list(mat.index), mat.values, list(mat.columns)
    genes = map(str, genes)
    mat = mat - mat.mean()
    with open(nf, "w") as f:
        column_header = "\t".join(['UNIQID', 'NAME', 'GWEIGHT'] + columns) + '\n'
        f.write(column_header)
        eweight = "\t".join( ['EWEIGHT', '', ''] + ['1'] * len(columns)) + '\n'  ### format column-flat-clusters for export
        f.write(eweight)
        for i, g in enumerate(genes):
            line = "\t".join([g] * 2 + ['1'] + map(str, mat[i])) + '\n'
            f.write(line)


for f in glob("*.txt"):
    print(f)
    txt2cdt(f)
