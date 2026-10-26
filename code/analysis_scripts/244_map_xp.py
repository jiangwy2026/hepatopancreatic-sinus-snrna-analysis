import pandas as pd, gzip, re, sys
sys.stdout.reconfigure(encoding="utf-8")
FEAT = r"D:\桌面\single cell\analysis_output\GCF_003789085.1_feature_table.txt.gz"
xps = {
 "JAG1":["XP_069984471"], "DLL1":["XP_069989473"], "NOTCH1":["XP_069992807"],
 "WNT5A":["XP_069977048","XP_069980464","XP_069985836","XP_069991869"],
 "FZD4":["XP_069971667","XP_069987369"], "FZD5":["XP_069969759"],
 "THBS4":["XP_027231548"], "LAMA1":["XP_069989416","XP_069989551"],
 "POSTN":["XP_027225409","XP_069998165"], "VWF":["XP_069988451"],
 "ITGAV":["XP_069969627","XP_069986831"],
 "ITGB1":["XP_027239322","XP_069973256","XP_069973257","XP_069983458"],
 "ITGA4":["XP_069974570","XP_069986011"]}
want = {x for v in xps.values() for x in v}
m = {}
with gzip.open(FEAT,"rt",encoding="utf-8",errors="ignore") as f:
    hdr=f.readline().lstrip("#").strip().split("\t")
    for line in f:
        p=line.rstrip("\n").split("\t")
        if len(p)<len(hdr): continue
        d=dict(zip(hdr,p))
        if d.get("feature")!="mRNA": continue
        acc=d.get("related_accession","").split(".")[0]
        if acc in want: m[acc]=(d.get("GeneID",""),d.get("name",""))
print("mapped",len(m),"of",len(want))
for g,v in xps.items():
    print(g,[(x,m.get(x,("UNMAPPED",""))) for x in v])
