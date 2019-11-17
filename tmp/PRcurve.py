#!/usr/bin/env python2.7
#gff file pos: 1-based close end
#gene2refseq file pos: 0-based close end
#all coordinate converted to 0-based open end format like bed
import matplotlib as mpl
import numpy as np
import os

mpl.use("pdf")
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["figure.figsize"] = (4, 2.75)
mpl.rcParams["figure.dpi"] = 100
mpl.rcParams["savefig.transparent"] = True
mpl.rcParams["savefig.bbox"] = "tight"
mpl.rcParams["font.size"] = 10.0
mpl.rcParams["font.sans-serif"] = "Arial"
mpl.rcParams["savefig.format"] = "pdf"

#choose seaborn 'Set3' palette
#colors = [(0.5529412031173706, 0.8274509906768799, 0.7803921699523926), (0.9964936564950382, 0.9986466744366814, 0.7025759500615737), (0.7490965163006502, 0.7336563019191518, 0.8525028922978569), (0.9786851223777322, 0.5073126014541177, 0.45665513592607837), (0.5077278263428631, 0.694256073587081, 0.8222376111675712), (0.9910188394434312, 0.7065282758544473, 0.3844213832827175)]

import brewer2mpl
colors =  brewer2mpl.get_map( 'Set2', 'qualitative',8 ).mpl_colors
colors.extend(brewer2mpl.get_map( 'Set3', 'qualitative',8 ).mpl_colors)

import sys, re, gzip, joblib
from glob import glob
import pandas as pd
#import seaborn as sns
import pylab


def parsePgl(fn):
    loops = []
    with open(fn, "r") as f:
        header = f.readline()
        for line in f.readlines():
            line = line.rstrip("\n").split("\t")
            lchr, lstart, lend, rchr, rstart, rend, loopid, score = line[0:8]
            lstart, lend, rstart, rend, score = int(lstart), int(lend), int(rstart), int(rend), float(score)
            if lstart + lend > rstart + rend:
                lstart, lend, rstart, rend = rstart, rend, lstart, lend
            loops.append([lchr, lstart, lend, rchr, rstart, rend, score])
    if loops:
        loops.sort(key=lambda x: x[6], reverse=False)
    return loops

def parseBedpe(fn):
    loops = {}
    with open(fn, "r") as f:
        for line in f.readlines():
            if re.match(r"^#", line):
                continue
            line = line.rstrip("\n").split("\t")
            lchr, lstart, lend, rchr, rstart, rend = line[0:6]
            lstart, lend, rstart, rend = int(lstart), int(lend), int(rstart), int(rend)
            if lstart + lend > rstart + rend:
                lstart, lend, rstart, rend = rstart, rend, lstart, lend
            if lchr not in loops:
                loops[lchr] = []
            loops[lchr].append([lchr, lstart, lend, rchr, rstart, rend, 0])
    return loops

def checkOverlap(loopA, loopB):
    if not loopA[1] >= loopB[2] and not loopA[2] <= loopB[1] and not loopA[4] >= loopB[5] and not loopA[5] <= loopB[4]:
        return True
    else:
        return False

for sample in ["GM12878", "K562"]:
    for comp_type in glob(os.path.join("loops", "*")):
        if not os.path.isdir(comp_type):
            continue
        comp_type = os.path.basename(comp_type)
        print comp_type
        result_lines = {}
        result_file = os.path.join("loops", "%s_%s.txt" % (sample, comp_type))
        jd_file = os.path.join("loops", "%s_%s.jd" % (sample, comp_type))
        if not os.path.isfile(result_file) or not os.path.isfile(jd_file):
            with open(result_file, 'w') as rf:
                rf.write("%s\tHuman marked loops Discovered\tHuman marked loops number\tSensitivity\tPredicted loops overlapped with human marked loops\tPredicted loops\tPrecision\n" % sample)
                for f in glob(os.path.join("loops", comp_type, "%s_*.pgl" % sample)):
                    key = re.sub('.pgl$', '', f)
                    key = os.path.basename(key)
                    key = re.sub("^%s_" % sample, '', key)
                    loops = parsePgl(f)
                    validset = parseBedpe(os.path.join("humanmark", "%s_humanMark.bedpe" % sample))
                    total_num = 0
                    for chr in validset:
                        total_num += len(validset[chr])
                    line = []
                    TP = 0
                    TP_PR = 0
                    for chr in validset:
                        for valid_loop in validset[chr]:
                            valid_loop[6] = 0
                    for i, loop in enumerate(loops):
                        FP_flag = True
                        chr = loop[0]
                        for valid_loop in validset[chr]:
                            if checkOverlap(loop, valid_loop):
                                FP_flag = False
                                if valid_loop[6] == 0:
                                    TP += 1
                                    valid_loop[6] = 1
                        if not FP_flag:
                            TP_PR += 1
                        else:
                            pass
                        line.append([float(TP) / total_num, float(TP) / (i + 1)])
                    if len(loops) == 0:
                        line.append([0, 0])
                        i = -1
                    rf.write("\t".join([key, str(TP), str(total_num), str(line[-1][0]), str(TP_PR), str(i + 1), str(line[-1][1])]))
                    rf.write("\n")
                    result_lines[key] = line
            joblib.dump(result_lines, jd_file)
        else:
            result_lines = joblib.load(jd_file)
        ds = {}
        fig, ax = pylab.subplots(figsize=(4*0.8, 2.75*0.8))
        for i, key in enumerate(sorted(result_lines.keys())):
            linedata = pd.DataFrame(result_lines[key])
            TPR = linedata.loc[:, 0].values
            PPV = linedata.loc[:, 1].values
            if (len(TPR) > 200):
                new_TPR = []
                new_PPV = []
                for j in np.arange(0, len(TPR), float(len(TPR))/200):
                    new_TPR.append(TPR[int(j)])
                    new_PPV.append(PPV[int(j)])
                TPR = new_TPR
                PPV = new_PPV
            #ax.plot(TPR, PPV, color=colors[i], marker=",", linewidth=1, label=key)
            ax.plot(TPR, PPV, marker=",", linewidth=1, label=key)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        legend = ax.legend(prop={'size':6}, fancybox=False)
        frame = legend.get_frame()
        frame.set_edgecolor('silver')
        ax.set_xlim([0, 1])
        ax.set_ylim([0, 1])
        ax.set_xlabel("Recall (sensitivity)")
        ax.set_ylabel("Precision (PPV)")
        pylab.savefig(os.path.join("plots", "%s_%s.pdf" % (sample, comp_type)))

