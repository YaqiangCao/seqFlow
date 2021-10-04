import json 
from glob import glob
from copy import deepcopy


meta = {
    "type":"bigwig",
    "url":"http://137.187.134.252/data/caoy7/",
    "name":"",
    "options":{
            "color":"black",
        },
    "metadata":{
            "sample":"",
            "assay":"",
        }
}


def pre(rot="1.ShuaiLiu_Hi-TrAC/1.20210909_KZ2320/3.mouseTestATAC",sample="mouse",assay="ATAC-seq"):
    fs = glob("*.bw")
    fs.sort()
    data = []
    for f in fs:
        m = deepcopy( meta )
        n = f.split(".bw")[0]
        m["name"] = n
        m["url"] = m["url"]+rot+"/"+f
        m["metadata"]["sample"] = sample
        m["metadata"]["assay"] = assay
        data.append( m )
    with open("20210909_KZ2320_atac.json","w") as fo:
        json.dump(data,fo)

pre()
