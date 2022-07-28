import gzip
import click
import pandas as pd


def getUniDNA(f,metaf="../5.statFilter/rawCellStat_filter.txt"):
    meta = pd.read_csv(metaf,index_col=0,sep="\t")
    cells = set(meta.index)
    cs = {}
    nf = f.split("/")[-1].replace(".gz","")
    with open(nf,"w") as fo:
        for i, line in enumerate(gzip.open(f, "rt")):
            if i % 10000 == 0:
                print("%s reads parsed from %s" % (i, f))
            line = line.split("\n")[0].split("\t")
            rid = line[3]
            cid = "_".join(rid.split("_")[:2])
            if cid not in cells:
                continue
            umi = rid.split("_")[2]
            if cid not in cs:
                cs[cid] = { }
            r = (line[0], int(line[1]), int(line[2]))
            if r not in cs[cid]:
                cs[cid][r] = set([umi])
                fo.write("\t".join(line)+"\n")
            else:
                #same location same umi, continue
                if umi not in cs[cid][r]:
                    cs[cid][r].add(umi)
                    fo.write("\t".join(line)+"\n")
    
@click.command()
@click.option("-name", required=True, help="Sample name/id for the data.")
def main(name):
    getUniDNA("../3.DNA/"+name + "_DNA.bed.gz")


if __name__ == '__main__':
    main()
