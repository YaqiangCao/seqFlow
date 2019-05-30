import pandas as pd

ks = [2,0,4,1,3]
f = "tissue_degs_K_G5.kgg"
ss = {}
for i,line in enumerate(open(f)):
    if i == 0:
        continue
    line = line.split("\n")[0].split("\t")
    line[1] = int( line[1] )
    if line[1] not in ss:
        ss[line[1]] = []
    ss[ line[1] ].append( line[0] )
gs = {}
ns = []
for i,k in enumerate(ks):
    ns.extend( ss[k] )
    for g in ss[k]:
        gs[g] = i
    ts = [t.split("|")[1] for t in ss[k]]
    with open("c%s.list"%i,"w") as fo:
        fo.write("\n".join(ts))
gs = pd.Series(gs)
gs = gs.sort_values()
gs.to_csv("clusters.txt",sep="\t")



ds = {}
lines = []
for i, line in enumerate(open("tissue_degs_K_G5.cdt")):
    if i < 2:
        lines.append( line )
    else:
        t = line.split("\n")[0].split("\t")[0]
        ds[t] = line
with open("degs_clusters.cdt","w") as fo:
    fo.write("".join(lines))
    for i in ns:
        fo.write(ds[i])
