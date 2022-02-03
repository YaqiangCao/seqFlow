import json 
from glob import glob
from copy import deepcopy

metaL = {
    "type":"longrange",
    "url":"http://137.187.134.252:8081/data/caoy7/",
    "name":"",
    "options":{
            "color":"black",
        },
    "metadata":{
            "sample":"",
            "assay":"",
        }
}


metaW = {
    "type":"bigwig",
    "url":"http://137.187.134.252:8081/data/caoy7/",
    "name":"",
    "options":{
            "color":"black",
        },
    "metadata":{
            "sample":"",
            "assay":"",
        }
}

def pre(rot="1.ShuaiLiu_Hi-TrAC/19.20220128_KZ2373",sample="human"):
    data = []

    fs = glob("*_PETs_washU.txt.gz")
    fs.sort()
    for f in fs:
        m = deepcopy( metaL )
        n = f.split("_PETs_washU.txt.gz")[0]
        m["name"] = n
        m["url"] = m["url"]+rot+"/"+f
        m["metadata"]["sample"] = sample
        m["metadata"]["assay"] = "Hi-TrAC"
        data.append( m )

    fs = glob("*.bw")
    fs.sort()
    for f in fs:
        m = deepcopy( metaW )
        n = f.split(".bw")[0]
        m["name"] = n
        m["url"] = m["url"]+rot+"/"+f
        m["metadata"]["sample"] = sample
        m["metadata"]["assay"] = "ATAC-seq"
        data.append( m )
    with open("20220128_KZ2373.json","w") as fo:
        json.dump(data,fo)

pre()
