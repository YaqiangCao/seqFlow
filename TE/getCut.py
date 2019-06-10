import pandas as pd
cut = 0.8
mat = pd.read_table( "mm10Rep.txt",index_col=0 )
s = mat[ "mappability_36" ]
s = s[ s>=cut ]
mat = mat.loc[ s.index,: ]
mat.to_csv( "mm10Rep_cut.txt",sep="\t" )
