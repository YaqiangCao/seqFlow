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


def pre(rot="6.jeff/1.20201020",sample="mouse",assay="ChIP"):
    fs = glob("*.bw")
    fs.sort()
    data = []
    for f in fs:
        m = deepcopy( meta )
        n = f.split(".")[0]
        m["name"] = n
        m["url"] = m["url"]+rot+"/"+f
        m["metadata"]["sample"] = sample
        m["metadata"]["assay"] = assay
        data.append( m )
    with open("JeffLab_20201020.json","w") as fo:
        json.dump(data,fo)

pre()
