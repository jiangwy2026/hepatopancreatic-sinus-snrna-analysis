import os, numpy as np, pandas as pd, anndata as ad, scanpy as sc
from scipy.stats import spearmanr
O=r"D:\桌面\single cell\analysis_output"
# expression metrics
A=ad.read_h5ad(os.path.join(O,"adata_C2_within.h5ad")); A.X=A.layers["counts"].copy(); sc.pp.normalize_total(A,target_sum=1e4); sc.pp.log1p(A)
gid2i={str(g):i for i,g in enumerate(A.var_names.astype(str))}; sym2gid={}
for gid,s in zip(A.var_names.astype(str),A.var["sym"].astype(str)): sym2gid.setdefault(s,gid)
g2m=A.obs.G2M_score.values.astype(float); pt=A.obs.dpt_pseudotime.values.astype(float); grp=A.obs.group.astype(str).values
cc=sym2gid.get('ccnd2'); ccv=A[:,gid2i[cc]].X.toarray().ravel()
early=(pd.qcut(A.obs.dpt_pseudotime.astype(float),3,labels=False).values==0)
ann=pd.read_csv(r"D:\桌面\single cell\all.genes.marker.stat.xls.xls",sep="\t",encoding="gb18030",low_memory=False)
fam=dict(zip(ann.Symbol.astype(str),ann.TF_family.astype(str)))
sel=['br','Trl','Cfp1','Taf7','Kdm4B','Scm','foxo','Stat5b','Su(var)205','pnt']
rows=[]
for tf in sel:
 gid=sym2gid.get(tf)
 if not gid: continue
 v=A[:,gid2i[gid]].X.toarray().ravel()
 rows.append({"TF":tf,"family":fam.get(tf,''),"pct_detected":100*(v>0).mean(),"mean_logexpr":v.mean(),
              "rho_G2M":spearmanr(v,g2m).correlation,"rho_ccnd2":spearmanr(v,ccv).correlation,"rho_pseudotime":spearmanr(v,pt).correlation,
              "early_blank":v[early&(grp=='blank')].mean(),"early_treated":v[early&(grp=='treated')].mean()})
E=pd.DataFrame(rows)
# merge target evidence
rk=pd.read_csv(os.path.join(O,"PH_TF_target_rank_enrichment.csv"))[["TF","mean_rank_pct","top10pct_frac"]]
gl=pd.read_csv(os.path.join(O,"PH_TF_target_G2M_high_low.csv"))[["TF","G2M_high_minus_low","p"]].rename(columns={"p":"target_gl_p"})
en=pd.read_csv(os.path.join(O,"PH_TF_target_enrichment_up.csv"))[["TF","overlap_up","hypergeom_p","neglog10p"]]
E=E.merge(rk,on="TF",how="left").merge(gl,on="TF",how="left").merge(en,on="TF",how="left")
E["early_logFC_proxy"]=np.log2((E.early_treated+0.05)/(E.early_blank+0.05))
E.to_csv(os.path.join(O,"PH_TF_selected_evidence.csv"),index=False,encoding="utf-8-sig")
print(E.round(3).to_string(index=False))
