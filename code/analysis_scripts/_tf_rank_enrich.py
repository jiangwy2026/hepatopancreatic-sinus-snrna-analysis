import os, numpy as np, pandas as pd
O=r"D:\桌面\single cell\analysis_output"
D=pd.read_csv(os.path.join(O,"PH_early_treated_vs_blank_DE.csv"))
D=D[np.isfinite(D.scores)].copy(); D["names"]=D.names.astype(str)
D=D.sort_values("scores",ascending=False).reset_index(drop=True)
D["rankpct"]=1-D.index/(len(D)-1)
rank=dict(zip(D.names,D.rankpct))
adj=pd.read_csv(os.path.join(O,"scenic_adjacencies.csv")); adj.columns=["TF","target","importance"]
ann=pd.read_csv(r"D:\桌面\single cell\all.genes.marker.stat.xls.xls",sep="\t",encoding="gb18030",low_memory=False)
sym2gid=dict(zip(ann.Symbol.astype(str),ann.GeneID.astype(str)))
adj["gid"]=adj.target.astype(str).map(sym2gid); adj=adj[adj.gid.isin(rank)].copy()
rows=[]
for tf,d in adj.groupby("TF"):
 d=d.sort_values("importance",ascending=False).drop_duplicates("gid").head(50)
 vals=[rank[g] for g in d.gid if g in rank]
 if len(vals)<10: continue
 rows.append({"TF":tf,"n":len(vals),"mean_rank_pct":float(np.mean(vals)),"median_rank_pct":float(np.median(vals)),"top10pct_frac":float(np.mean(np.array(vals)>=0.9))})
R=pd.DataFrame(rows).sort_values("mean_rank_pct",ascending=False)
print(R.head(40).round(3).to_string(index=False))
R.to_csv(os.path.join(O,"PH_TF_target_rank_enrichment.csv"),index=False,encoding="utf-8-sig")
