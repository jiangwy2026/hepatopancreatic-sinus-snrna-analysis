import os, numpy as np, pandas as pd, anndata as ad, scanpy as sc
from scipy.stats import mannwhitneyu
O=r"D:\桌面\single cell\analysis_output"
A=ad.read_h5ad(os.path.join(O,"adata_C2_within.h5ad")); A.X=A.layers["counts"].copy(); sc.pp.normalize_total(A,target_sum=1e4); sc.pp.log1p(A)
x=np.asarray(A.X.toarray(),dtype=np.float32); Xz=(x-x.mean(0,keepdims=True))/(x.std(0,keepdims=True)+1e-6)
gid2i={str(g):i for i,g in enumerate(A.var_names.astype(str))}
# symbol to gid
sym2gid={}
for gid,s in zip(A.var_names.astype(str),A.var["sym"].astype(str)): sym2gid.setdefault(s,gid)
adj=pd.read_csv(os.path.join(O,"scenic_adjacencies.csv")); adj.columns=["TF","target","importance"]
adj["gid"]=adj.target.astype(str).map(sym2gid); adj=adj[adj.gid.notna()].copy()
g2m=A.obs.G2M_score.values.astype(float)
hi=g2m>=np.quantile(g2m,0.75); lo=g2m<=np.quantile(g2m,0.25)
rows=[]
for tf,d in adj.groupby("TF"):
 d=d.sort_values("importance",ascending=False).drop_duplicates("gid").head(50)
 idx=[gid2i[g] for g in d.gid if g in gid2i]
 if len(idx)<10: continue
 s=Xz[:,idx].mean(1)
 u,p=mannwhitneyu(s[hi],s[lo],alternative="two-sided")
 diff=float(s[hi].mean()-s[lo].mean())
 rows.append({"TF":tf,"n_targets":len(idx),"G2M_high_minus_low":diff,"p":p})
R=pd.DataFrame(rows).sort_values("G2M_high_minus_low",ascending=False)
print(R.head(40).round(4).to_string(index=False))
R.to_csv(os.path.join(O,"PH_TF_target_G2M_high_low.csv"),index=False,encoding="utf-8-sig")
