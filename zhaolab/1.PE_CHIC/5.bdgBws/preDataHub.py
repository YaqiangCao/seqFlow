import json 
from glob import glob
from copy import deepcopy

metaL = {
    "type":"longrange",
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


metaW = {
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

def pre(rot="2.GuangzheGe_scDNase/2.20211021_KZ2232/2.mouse/",sample="mouse"):
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
        m["metadata"]["assay"] = "Hi-TrAC 1D"
        data.append( m )
    with open("20211021_KZ2232_mouse_R1.json","w") as fo:
        json.dump(data,fo)

pre()
