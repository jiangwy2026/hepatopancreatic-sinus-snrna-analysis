# -*- coding: utf-8 -*-
"""283_ph_marker_orthologs.py — PH marker shrimp-human BLAST + specificity."""
import os, gzip, subprocess, sys
import numpy as np, pandas as pd

OUT=r"D:\桌面\single cell\analysis_output"
BLAST=r"D:\桌面\single cell\tools\ncbi-blast-2.17.0+\bin"
W=os.path.join(os.environ["TEMP"],"phmarker"); os.makedirs(W,exist_ok=True)
markers={"PPAF2":["113812814","113823945"],"IGFBP7":["113811702"],"GRIK2":["113823802"],
         "SSPO/hemocytin":["113824171"],"CCND2":["113814652","113821812"],
         "TGM":["113823945","113823934","113823455"],"SLO":["113807409"],
         "SERPINB1":["113799939"],"SYNE1":["113823000"]}
# map GeneID -> XP
xp={}
with gzip.open(os.path.join(OUT,"GCF_003789085.1_feature_table.txt.gz"),"rt",encoding="utf-8",errors="ignore") as f:
    hdr=f.readline().lstrip("#").strip().split("\t")
    for line in f:
        p=line.rstrip("\n").split("\t")
        if len(p)<len(hdr): continue
        d=dict(zip(hdr,p))
        if d.get("feature")!="mRNA": continue
        xp[d.get("GeneID","")]=d.get("related_accession","").split(".")[0]
# extract shrimp seqs
want={}
for name,gids in markers.items():
    for g in gids:
        if g in xp: want[xp[g]]=name+"|"+g
seqs={}; keep=False; acc=None
for line in open(os.path.join(OUT,"GCF_003789085.1_protein.faa"),encoding="utf-8",errors="ignore"):
    if line.startswith(">"):
        acc=line[1:].split()[0].split(".")[0]; keep=acc in want
        if keep: seqs[acc]=[line]
    elif keep: seqs[acc].append(line)
with open(os.path.join(W,"shrimp_markers.faa"),"w") as f:
    for a,lines in seqs.items(): f.write(lines[0] if lines[0].startswith(">") else ">"+a+"\n"); f.writelines(lines[1:] if lines[0].startswith(">") else lines)
print("shrimp marker seqs",len(seqs))
# human db
hfaa=os.path.join(W,"human.faa")
if not os.path.exists(hfaa):
    with gzip.open(os.path.join(OUT,"human_refseq_protein.faa.gz"),"rt",encoding="utf-8",errors="ignore") as f, open(hfaa,"w") as o: o.write(f.read())
if not os.path.exists(os.path.join(W,"human_db.phr")):
    subprocess.run([os.path.join(BLAST,"makeblastdb.exe"),"-in",hfaa,"-dbtype","prot","-out",os.path.join(W,"human_db")],check=True,stdout=subprocess.DEVNULL)
subprocess.run([os.path.join(BLAST,"blastp.exe"),"-query",os.path.join(W,"shrimp_markers.faa"),"-db",os.path.join(W,"human_db"),
                "-evalue","1e-3","-max_target_seqs","3","-outfmt","6 qseqid sseqid pident length qlen slen evalue bitscore qcovs",
                "-out",os.path.join(W,"ph_human.tsv")],check=True,stdout=subprocess.DEVNULL)
b=pd.read_csv(os.path.join(W,"ph_human.tsv"),sep="\t",header=None,names=["q","s","pident","length","qlen","slen","evalue","bitscore","qcovs"])
b["name"]=b.q.map(lambda x: want.get(x.split(".")[0],"?"))
best=b.sort_values("bitscore",ascending=False).groupby("name").first().reset_index()
# human descriptions
desc={}
with gzip.open(os.path.join(OUT,"human_refseq_protein.faa.gz"),"rt",encoding="utf-8",errors="ignore") as f:
    for l in f:
        if l.startswith(">"):
            acc=l[1:].split()[0].split(".")[0]; desc[acc]=l[1:].strip()[len(l[1:].split()[0]):].strip()[:90]
best["human_desc"]=best.s.map(lambda x: desc.get(x.split(".")[0],""))
# specificity
mk=pd.read_csv(os.path.join(OUT,"markers_10groups.csv"))
best["marker"]=best.name.str.split("|").str[0]; best["gid"]=best.name.str.split("|").str[1]
rows=[]
for r in best.itertuples():
    hit=mk[(mk.grp=="PH")&(mk.names.astype(str)=="ncbi_"+r.gid)]
    if len(hit):
        h=hit.sort_values("scores",ascending=False).iloc[0]
        rows.append({"name":r.name,"marker":r.marker,"marker_gid":"ncbi_"+r.gid,"log2FC":h.logfoldchanges,"pct_PH":h.pct_nz_group*100,"pct_other":h.pct_nz_reference*100})
    else: rows.append({"name":r.name,"marker":r.marker,"marker_gid":"ncbi_"+r.gid})
sp=pd.DataFrame(rows)
out=best.merge(sp,on="name",how="left")
out.to_csv(os.path.join(OUT,"step1","PH_marker_human_orthologs.csv"),index=False,encoding="utf-8-sig")
pd.set_option("display.width",260); pd.set_option("display.max_colwidth",50)
print(out[["name","s","human_desc","pident","qcovs","evalue","log2FC","pct_PH","pct_other"]].to_string(index=False))
