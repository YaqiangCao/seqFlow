import gzip
import click
import pysam
import pandas as pd


def getUniRNA(f,name,metaf="../5.statFilter/rawCellStat_filter.txt"):
    meta = pd.read_csv(metaf,index_col=0,sep="\t")
    cells = set(meta.index)
    cs = {}
    nf = name + "_RNA.bam"
    with pysam.AlignmentFile(nf,"w",template=pysam.AlignmentFile(f)) as fo:
        for i,b in enumerate( pysam.AlignmentFile(f).fetch(until_eof=True) ) :
            if i % 10000 == 0:
                print("%s reads parsed from %s" % (i, f))
            rid = b.query_name
            cid = "_".join(rid.split("_")[:2])
            if cid not in cells:
                continue
            umi = rid.split("_")[2]
            if cid not in cs:
                cs[cid] = { }
            r = (b.reference_name, b.reference_start, b.reference_end)
            if r not in cs[cid]:
                cs[cid][r] = set([umi])
                fo.write(b)
            else:
                #same location same umi, continue
                if umi not in cs[cid][r]:
                    cs[cid][r].add(umi)
                    fo.write(b)
                
@click.command()
@click.option("-name", required=True, help="Sample name/id for the data.")
def main(name):
    getUniRNA("../4.RNA/%s_Aligned.out.bam"%name,name)


if __name__ == '__main__':
    main()
