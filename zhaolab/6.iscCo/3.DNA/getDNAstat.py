
import gzip
import click
import pandas as pd



 

def getStat(f):
    cs = {}
    for i,line in enumerate(gzip.open(f,"rt")):
        if i % 10000 == 0:
            print("%s reads parsed from %s"%(i,f))
        line = line.split("\n")[0].split("\t")
        rid = line[3]
        cid = "_".join( rid.split("_")[:2] )
        umi = rid.split("_")[2]
        if cid not in cs:
            cs[cid] = {
                    "total": 0,
                    "reads": {},
                }
        cs[cid]["total"] += 1
        r = (line[0],int(line[1]),int(line[2]))
        if r not in cs[cid]["reads"]:
           cs[cid]["reads"][r] = set([umi])
        else:
            #same location same umi, continue
            if umi in cs[cid]["reads"][r]:
                continue
            else:
               cs[cid]["reads"][r].add(umi) 
    ds = {}
    for k,v in cs.items():
        c = 0
        #count read and umi
        for r, u in v["reads"].items():
            c += len(u)
        ds[ k ] = {
            "totalMapped(mapq>=10)":v["total"],
            "uniqueMapped": c,
            "redundancy": 1 - c/v["total"]
            }
    ds = pd.DataFrame(ds).T
    s = list(ds.index)
    s.sort()
    ds = ds.loc[s,]
    ds.to_csv("cellStat.txt",sep="\t",index_label="cellId")



@click.command()
@click.option("-name",
              required=True,
              help="Sample name/id for the data.")
def main(name):
    getStat(name+"_DNA_all.bed.gz")

if __name__ == '__main__':
    main()
