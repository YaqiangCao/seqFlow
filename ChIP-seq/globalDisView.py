#--coding:utf-8--
"""
globalDisView.py global data distribution view
"""

__author__ = "CAO Yaqiang"
__date__ = "2019-05-17"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#sys
from glob import glob
from datetime import datetime
import warnings
warnings.filterwarnings("ignore")

#3rd
import matplotlib as mpl
mpl.use("pdf")
import seaborn as sns
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["figure.figsize"] = (4, 2.75)
mpl.rcParams["figure.dpi"] = 100
mpl.rcParams["savefig.transparent"] = True
mpl.rcParams["savefig.bbox"] = "tight"
mpl.rcParams["font.size"] = 10.0
mpl.rcParams["font.sans-serif"] = "Arial"
mpl.rcParams["savefig.format"] = "pdf"
import pylab
sns.set_style("white")
import numpy as np
import pandas as pd
import brewer2mpl
colors = brewer2mpl.get_map('Set1', 'qualitative', 9).mpl_colors
colors.extend(brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors)
del colors[5]
del colors[10]


def kdePlot(data, pre, name):
    """
    KDE plot for checking 
    """
    #density plot
    f, ax = pylab.subplots()
    for c in data.columns:
        d = data[c]
        sns.kdeplot(d, label=c, shade=False)
    leg = ax.legend(loc="best", shadow=True, fancybox=True)
    pylab.setp(leg.get_texts())
    ax.set_xlabel(name)
    ax.set_ylabel("Density")
    #ax.set_xlim((15,30))
    pylab.tight_layout()
    pylab.savefig(pre + "_kde.pdf", bbox_inches='tight')


def rlePlot(data, pre, name):
    """
    Relative log expression (RLE) plots.
    Reference to:https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0191629 
    """
    print(data.shape)
    data = data.sub(data.median(axis=1), axis='index')
    fig, ax = pylab.subplots(figsize=(4, 2.75))
    sns.boxplot(data=data, fliersize=0, color=colors[0])
    ax.set_ylim([-2, 2])
    ax.set_ylabel("RLE")
    ax.set_xticklabels(list(data.columns), rotation=45, fontsize=8, ha="right")
    ax.set_title(pre)
    pylab.tight_layout()
    pylab.savefig(pre + "_kde.pdf", bbox_inches='tight')


def displot():
    fs = glob("../10.getNoiseFilter/*.txt")
    for f in fs:
        pre = f.split("/")[-1].split(".")[0]
        data = pd.read_table(f, index_col=0)
        data = np.log2(data + 1.0)
        kdePlot(data, pre, "log2(TPM)")
        rlePlot(data, pre, "log2(TPM)")


def main():
    displot()


main()
