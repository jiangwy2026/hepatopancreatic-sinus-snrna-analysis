import os, numpy as np, pandas as pd
from scipy.stats import hypergeom
O=r"D:\桌面\single cell\analysis_output"
D=pd.read_csv(os.path.join(O,"PH_early_treated_vs_blank_DE.csv"))
D=D[np.isfinite(D.logfoldchanges)].copy(); D["names"]=D.names.astype(str)
up=D[(D.pvals_adj<0.1)&(D.logfoldchanges>0)].sort_values("logfoldchanges",ascending=False).head(500).names.tolist()
if len(up)<50: up=D.sort_values("logfoldchanges",ascending=False).head(500).names.tolist()
background=set(D.names)
adj=pd.read_csv(os.path.join(O,"scenic_adjacencies.csv")); adj.columns=["TF","target","importance"]
ann=pd.read_csv(r"D:\桌面\single cell\all.genes.marker.stat.xls.xls",sep="\t",encoding="gb18030",low_memory=False)
sym2gid=dict(zip(ann.Symbol.astype(str),ann.GeneID.astype(str)))
adj["gid"]=adj.target.astype(str).map(sym2gid)
adj=adj[adj.gid.isin(background)].copy()
print("adj mapped",len(adj),"up",len(up))
rows=[]
for tf,d in adj.groupby("TF"):
    d=d.sort_values("importance",ascending=False).drop_duplicates("gid").head(50)
    targets=set(d.gid); N=len(targets); n=len(up); k=len(targets & set(up))
    if N<10: continue
    p=hypergeom.sf(k-1,len(background),N,n)
    rows.append({"TF":tf,"n_targets":N,"overlap_up":k,"hypergeom_p":p,"neglog10p":-np.log10(max(p,1e-300))})
R=pd.DataFrame(rows).sort_values(["neglog10p","overlap_up"],ascending=False)
print(R.head(40).round(3).to_string(index=False))
R.to_csv(os.path.join(O,"PH_TF_target_enrichment_up.csv"),index=False,encoding="utf-8-sig")
