import gzip 
from glob import glob
from joblib import Parallel,delayed

def txt2bedpe(f):
    print f
    fn = f.split("/")[-1].replace(".txt",".bedpe")
    with open(fn,"w") as fo:
        #for i,line in enumerate(gzip.open(f)):
        for i,line in enumerate(open(f)):
            line = line.split( "\n" )[ 0 ].split( "\t" ) 
            nline = line[0:6]
            nline.extend([i,".",line[7],line[8]])
            fo.write("\t".join(map(str,nline))+"\n")

fs = glob("*.txt")
Parallel( n_jobs=len(fs) )( delayed( txt2bedpe )( f ) for f in fs )
