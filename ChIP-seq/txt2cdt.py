import glob
import pandas as pd
import numpy as np
import string


def txt2cdt(f):
    mat = pd.read_table(f, index_col=0)
    genes, mat, columns = list(mat.index), mat.values, list(mat.columns)
    genes = map(str, genes)
    mat = mat - mat.mean()
    nf = f.replace(".txt", ".cdt")
    with open(nf, "w") as f:
        column_header = string.join(['UNIQID', 'NAME', 'GWEIGHT'] + columns,
                                    '\t') + '\n'
        f.write(column_header)
        eweight = string.join(
            ['EWEIGHT', '', ''] + ['1'] * len(columns),
            '\t') + '\n'  ### format column-flat-clusters for export
        f.write(eweight)
        for i, g in enumerate(genes):
            line = string.join([g] * 2 + ['1'] + map(str, mat[i]), '\t') + '\n'
            f.write(line)


for f in glob.glob("*.txt"):
    print f
    txt2cdt(f)
