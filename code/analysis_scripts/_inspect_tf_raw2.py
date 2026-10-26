import anndata as ad, numpy as np, pandas as pd, os, scanpy as sc
from scipy.stats import spearmanr
O=r"D:\桌面\single cell\analysis_output"
A=ad.read_h5ad(os.path.join(O,"adata_C2_within.h5ad")); A.X=A.layers["counts"].copy(); sc.pp.normalize_total(A,target_sum=1e4); sc.pp.log1p(A)
gid2i={str(g):i for i,g in enumerate(A.var_names.astype(str))}; sym2gid={}
for gid,s in zip(A.var_names.astype(str),A.var["sym"].astype(str)): sym2gid.setdefault(s,gid)
g2m=A.obs.G2M_score.values.astype(float); pt=A.obs.dpt_pseudotime.values.astype(float); grp=A.obs.group.astype(str).values; early=(pd.qcut(A.obs.dpt_pseudotime.astype(float),3,labels=False).values==0)
cc=sym2gid.get('ccnd2'); y=A[:,gid2i[cc]].X.toarray().ravel()
for tf in ['pnr','tim','br','Trl','Cfp1','Taf7','Kdm4B','pnt','Stat5b']:
 gid=sym2gid.get(tf)
 if not gid: print(tf,'not found'); continue
 v=A[:,gid2i[gid]].X.toarray().ravel()
 print(tf,'pct',round((v>0).mean()*100,1),'mean',round(v.mean(),3),'rhoG2M',round(spearmanr(v,g2m).correlation,3),'rhoCC',round(spearmanr(v,y).correlation,3),'rhoPT',round(spearmanr(v,pt).correlation,3),'earlyB/T',round(v[early&(grp=='blank')].mean(),2),round(v[early&(grp=='treated')].mean(),2))
