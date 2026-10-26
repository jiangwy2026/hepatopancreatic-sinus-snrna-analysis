# -*- coding: utf-8 -*-
"""199_focused_lr_screen.py - candidate orthology-based ligand-receptor screen."""
import os, sys
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
OUT=r"D:\桌面\single cell\analysis_output"
GROUPS={"0":"TS","1":"GH","2":"PH","3":"MH2","4":"PM-like","5":"F","6":"R","7":"B","8":"MyoEpi","9":"M"}
ORDER=["PH","GH","TS","MH2","F","R","B","M","MyoEpi","PM-like"]
SENDERS=["F","R","B","M","MyoEpi","PM-like"]; RECEIVERS=["PH","GH","TS","MH2"]
PAIRS=[
 ("Wnt5b","ncbi_113806005","Fzd4","ncbi_113813183","Wnt"),
 ("Wnt6","ncbi_113812192","Fzd1","ncbi_113813182","Wnt"),
 ("PDGF/VEGF-like","ncbi_113804544","Flt1","ncbi_113803634","VEGF/PDGF"),
 ("FGF16","ncbi_113829301","FGFR1","ncbi_113803204","FGF"),
 ("THBS4","ncbi_113824160","ITGA5","ncbi_113807444","ECM-integrin"),
 ("THBS4","ncbi_113824160","ITGB3","ncbi_113805718","ECM-integrin"),
 ("THBS4","ncbi_113824160","SDC2","ncbi_113807065","ECM-syndecan"),
 ("Igfbp7","ncbi_113811702","Igf1r","ncbi_113817108","IGF"),
]
def main():
 A=ad.read_h5ad(os.path.join(OUT,"adata_soupx.h5ad")); A.X=A.layers["counts"].copy(); sc.pp.normalize_total(A,target_sum=1e4); sc.pp.log1p(A)
 ct=A.obs.vendor_cluster.astype(str).map(GROUPS).values; gid2i={str(g):i for i,g in enumerate(A.var_names.astype(str))}
 rows=[]
 for lname,lgid,rname,rgid,path in PAIRS:
  if lgid not in gid2i or rgid not in gid2i: continue
  lv=A[:,gid2i[lgid]].X.toarray().ravel(); rv=A[:,gid2i[rgid]].X.toarray().ravel()
  for s in SENDERS:
   ms=ct==s; lm=float(lv[ms].mean()); lp=float((lv[ms]>0).mean()*100)
   for r in RECEIVERS:
    mr=ct==r; rm=float(rv[mr].mean()); rp=float((rv[mr]>0).mean()*100)
    if lp<5 or rp<5: continue
    rows.append({"pathway":path,"ligand":lname,"receptor":rname,"sender":s,"receiver":r,"ligand_mean":lm,"ligand_pct":lp,"receptor_mean":rm,"receptor_pct":rp,"score":np.sqrt(lm*rm)})
 R=pd.DataFrame(rows).sort_values("score",ascending=False)
 R.to_csv(os.path.join(OUT,"focused_LR_candidates.csv"),index=False,encoding="utf-8-sig")
 print(R.head(30).round(3).to_string(index=False))
 print("pairs",len(R))
if __name__=="__main__": main()
