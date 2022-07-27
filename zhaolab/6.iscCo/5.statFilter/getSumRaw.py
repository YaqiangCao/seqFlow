from collections import Counter
import click
import numpy as np
import pandas as pd
from pipi.sets import *


def get(fa="../2.sepSingleDNARNA/cellStat.txt",
        fb="../3.DNA/cellStat.txt",
        fc="../4.RNA/cellStat.txt"):
    mata = pd.read_csv(fa, index_col=0, sep="\t")
    matb = pd.read_csv(fb, index_col=0, sep="\t")
    matc = pd.read_csv(fc, index_col=0, sep="\t")
    mat = {}
    mat["1_totalDNA"] = mata["totalDNA"]
    mat["2_totalRNA"] = mata["totalRNA"]
    mat["3_mappedDNA(mapq>=10)"] = matb["totalMapped(mapq>=10)"]
    mat["4_uniqueDNA"] = matb["uniqueMapped"]
    mat["5_DNARedundancy"] = matb["redundancy"]
    mat["6_mappedRNA"] = matc["totalMapped"]
    mat["7_uniqueRNA"] = matc["uniqueMapped"]
    mat["8_RNARedundancy"] = matc["redundancy"]
    mat = pd.DataFrame(mat)
    mat = mat.fillna(0.0)
    mat.to_csv("rawCellStat.txt", sep="\t", index_label="cellId")


def plot():
    mat = pd.read_csv("rawCellStat.txt", sep="\t", index_col=0)
    sa = mat["4_uniqueDNA"]
    sa = sa[sa > 0]
    sb = mat["7_uniqueRNA"]
    sb = sb[sb > 0]
    s = sa.index.intersection(sb.index)
    sa = np.log10(sa[s])
    sb = np.log10(sb[s])
    fig, ax = pylab.subplots()
    data = pd.DataFrame({"DNA reads, log10": sa, "RNA reads, log10": sb})
    #g = sns.jointplot(x=sa,y=sb,kind="kde",fill=True,marginal_ticks=True)
    g = sns.jointplot(data=data,
                      x="DNA reads, log10",
                      y="RNA reads, log10",
                      kind="kde",
                      fill=True,
                      marginal_ticks=True)
    pylab.savefig("rawCellStat.pdf")


def filterCell(dna=1000, rna=100):
    mat = pd.read_csv("rawCellStat.txt", sep="\t", index_col=0)
    sa = mat["4_uniqueDNA"]
    sa = sa[sa > dna]
    sb = mat["7_uniqueRNA"]
    sb = sb[sb > rna]
    s = sa.index.intersection(sb.index)
    mat = mat.loc[s, ]
    mat.to_csv("rawCellStat_filter.txt", sep="\t", index_label="cellId")
    cs = {}
    for t in s:
        nc = t.split("_")
        w = nc[0]
        c = nc[1]
        if w not in cs:
            cs[w] = set()
        cs[w].add(c)
    for k, v in cs.items():
        cs[k] = len(v)
    cs = pd.Series(cs)
    fig, ax = pylab.subplots()
    sns.histplot(x=cs, binwidth=1, ax=ax)
    ax.set_xlabel(
        "number of cells in each well\n DNA reads>%s, RNA reads >%s" %
        (dna, rna))
    pylab.savefig("numberOfCellsPerWell.pdf")


def getSum(name,
           fa="../2.sepSingleDNARNA/summary.txt",
           fb="rawCellStat.txt",
           fc="rawCellStat_filter.txt"):
    mata = pd.read_csv(fa, index_col=0, sep="\t")
    sa = mata.sum()
    matb = pd.read_csv(fb, index_col=0, sep="\t")
    sb = matb.sum()
    matc = pd.read_csv(fc, index_col=0, sep="\t")
    sc = matc.sum()
    rs = {
        "1_totalRawReads": sa["total"],
        "2_withBarcodeReads": sa["withBarcode"],
        "3_withLinkerReads": sa["withLinker"],
        "4_seperatedDNARawReads": sa["DNA"],
        "5_seperatedRNARawReads": sa["RNA"],
        "6_mappedDNA(mapq>=10)": sb["3_mappedDNA(mapq>=10)"],
        "7_DNAMappingRatio": sb["3_mappedDNA(mapq>=10)"] / sa["DNA"],
        "8_DNAUniqueReads": sb["4_uniqueDNA"],
        "9_DNAYield": sb["4_uniqueDNA"] / sa["DNA"],
        "10_filteredDNAReads": sc["4_uniqueDNA"],
        "11_filteredDNAYield": sc["4_uniqueDNA"] / sa["DNA"],
        "12_mappedRNA": sb["6_mappedRNA"],
        "13_RNAMappingRatio": sb["6_mappedRNA"] / sa["RNA"],
        "14_RNAYield": sb["7_uniqueRNA"] / sa["RNA"],
        "15_filteredRNAReads": sc["7_uniqueRNA"],
        "16_filteredRNAYield": sc["7_uniqueRNA"] / sa["RNA"],
        "17_medianDNARedundancy": np.median(matc["5_DNARedundancy"]),
        "18_medianRNARedundancy": np.median(matc["8_RNARedundancy"]),
    }
    rs = pd.Series(rs)
    rs.to_csv(name + "_statReport.txt", sep="\t",header=None)


@click.command()
@click.option("-name", required=True, help="Sample name/id for the data.")
def main(name):
    get()
    plot()
    filterCell()
    getSum(name)


main()
