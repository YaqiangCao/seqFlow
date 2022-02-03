import os
from glob import glob
from joblib import Parallel, delayed


def get(f):
    n = f.split("/")[-1].split(".bam")[0]
    if os.path.isfile(n+".bw"):
        return 
    cmd = "bamCoverage -b %s -o %s.bw --ignoreDuplicates --minMappingQuality 10 --normalizeUsing CPM"%(f, n)
    print(cmd)
    os.system(cmd)

fs = glob("../2.mapping/*/*.bam")
Parallel(n_jobs=min(len(fs),20))(delayed(get)(f) for f in fs)
