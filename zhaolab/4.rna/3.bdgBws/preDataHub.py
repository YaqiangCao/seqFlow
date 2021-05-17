import json 
from glob import glob
from copy import deepcopy

meta = {
    "type":"bigwig",
    "url":"https://hpc.nih.gov/~caoy7/",
    "name":"",
    "options":{
            "color":"black",
        },
    "metadata":{
            "sample":"",
            "assay":"",
        }
}


def pre(rot="1.Trac/tmp/8.20210514_KZ2228/2.k562_4tfs_rna",sample="mouse",assay="RNA-seq"):
    fs = glob("*.bw")
    fs.sort()
    data = []
    for f in fs:
        m = deepcopy( meta )
        n = f.split(".")[0]
        m["name"] = n
        m["url"] = m["url"]+rot+"/"+f
        m["metadata"]["sample"] = sample
        m["metadata"]["assay"] = n.split("_")[0]
        data.append( m )
    with open("20210514_KZ2228_K562_RNA.json","w") as fo:
        json.dump(data,fo)

pre()
