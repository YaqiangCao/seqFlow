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


def pre(rot="13.GuangZhe_naiveCD4/1.20210422_KZ2218",sample="human",assay="iscDNase-seq"):
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
    with open("20210422_KZ2218.json","w") as fo:
        json.dump(data,fo)

pre()
