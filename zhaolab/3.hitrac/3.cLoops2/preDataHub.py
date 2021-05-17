import json 
from glob import glob
from copy import deepcopy

meta = {
    "type":"longrange",
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


def pre(rot="1.Trac/tmp/6.KZ2214/1.human_hitrac",sample="mouse",assay="Hi-TrAC"):
    fs = glob("*_PETs_washU.txt.gz")
    fs.sort()
    data = []
    for f in fs:
        m = deepcopy( meta )
        n = f.split("_PETs_washU.txt.gz")[0]
        m["name"] = n
        m["url"] = m["url"]+rot+"/"+f
        m["metadata"]["sample"] = sample
        m["metadata"]["assay"] = assay
        data.append( m )
    with open("KZ2214_human_hitrac.json","w") as fo:
        json.dump(data,fo)

pre()
