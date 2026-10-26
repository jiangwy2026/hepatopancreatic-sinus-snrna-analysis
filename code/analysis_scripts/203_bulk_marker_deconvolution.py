# -*- coding: utf-8 -*-
"""203_bulk_marker_deconvolution.py - marker-based estimate from full bulk matrix."""
import os, sys
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
from scipy.optimize import nnls
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
OUT=r"D:\桌面\single cell\analysis_output"; PREV=os.path.join(OUT,"figures_pub","sc_summary_new_names","individual_previews")
GROUPS={"0":"TS","1":"GH","2":"PH","3":"MH2","4":"PM-like","5":"F","6":"R","7":"B","8":"MyoEpi","9":"M"}
ORDER=["PH","GH","TS","MH2","F","R","B","M","MyoEpi","PM-like"]
def main():
 B=pd.read_csv(os.path.join(OUT,"bulk_GHRL25030101","5.ExpressionStat","all.genes.expression.annot.xls"),sep="\t",usecols=["Symbol","haemospasia_count","untreated_count"])
 B["Symbol"]=B.Symbol.astype(str); B=B[B.Symbol!="--"].copy(); B=B.groupby("Symbol",as_index=False).agg(untreated_count=("untreated_count","sum"),haemospasia_count=("haemospasia_count","sum")).set_index("Symbol")
 A=ad.read_h5ad(os.path.join(OUT,"adata_soupx.h5ad")); A.X=A.layers["counts"].copy(); sc.pp.normalize_total(A,target_sum=1e4); sc.pp.log1p(A)
 ct=A.obs.vendor_cluster.astype(str).map(GROUPS).values; sym2gid={}
 for gid,s in zip(A.var_names.astype(str),A.var["sym"].astype(str)): sym2gid.setdefault(s,gid)
 M=pd.read_csv(os.path.join(OUT,"top_markers_10groups.csv")); M=M[M["rank"]<=10]
 rows=[]
 for g in ORDER:
  syms=sorted(set(M[M.grp==g].sym.astype(str)) & set(B.index))
  d=B.loc[syms] if syms else B.iloc[0:0]
  if len(syms)>=1:
   s0=float(np.log1p(d.untreated_count).mean()); s6=float(np.log1p(d.haemospasia_count).mean())
  else: s0=s6=np.nan
  rows.append({"group":g,"n_markers":len(syms),"marker_score_untreated":s0,"marker_score_haemospasia":s6})
 R=pd.DataFrame(rows)
 for c in ["marker_score_untreated","marker_score_haemospasia"]:
  v=R[c].clip(lower=0); R[c+"_norm"]=v/(v.sum()+1e-9)
 # marker-only NNLS on union of markers
 syms=sorted(set(M.sym.astype(str)) & set(B.index))
 gids=[sym2gid[s] for s in syms if s in sym2gid]
 # ensure unique symbols map to unique gids
 smap={s:sym2gid[s] for s in syms if s in sym2gid}
 syms=list(smap); X=np.asarray(A[:,[smap[s] for s in syms]].X.toarray(),dtype=float)
 S=np.vstack([X[ct==g].mean(0) for g in ORDER]).T; S=S/(S.sum(0,keepdims=True)+1e-9)
 b0=B.loc[syms,"untreated_count"].values.astype(float); b6=B.loc[syms,"haemospasia_count"].values.astype(float)
 b0=b0/(b0.sum()+1e-9); b6=b6/(b6.sum()+1e-9)
 def fit(b):
  c,_=nnls(S,b); return c/(c.sum()+1e-9)
 R["nnls_marker_untreated"]=fit(b0); R["nnls_marker_haemospasia"]=fit(b6)
 R.to_csv(os.path.join(OUT,"bulk_marker_deconvolution_results.csv"),index=False,encoding="utf-8-sig")
 print(R.round(4).to_string(index=False)); print("markers",len(syms))
 fig=plt.figure(figsize=(5.2,3.4),facecolor="white"); ax=fig.add_axes([0.20,0.20,0.77,0.70]); y=np.arange(len(R))
 ax.barh(y-0.18,R.marker_score_untreated_norm,height=0.34,color="#4C78A8",label="Untreated")
 ax.barh(y+0.18,R.marker_score_haemospasia_norm,height=0.34,color="#D6604D",label="Haemospasia")
 ax.set_yticks(y); ax.set_yticklabels(R.group,fontsize=5.6); ax.invert_yaxis(); ax.set_xlabel("Marker-based proportion",fontsize=5.6)
 ax.legend(frameon=False,fontsize=5.0,loc="lower right"); ax.set_title("Bulk marker deconvolution",fontsize=6.8,pad=4)
 for sp in ("top","right"): ax.spines[sp].set_visible(False)
 ax.tick_params(axis="x",labelsize=4.8,length=1.0); ax.tick_params(axis="y",length=0)
 p=os.path.join(PREV,"06c_bulk_marker_deconvolution_preview.png"); fig.savefig(p,dpi=400,bbox_inches="tight",pad_inches=0.04,facecolor="white"); plt.close(fig); print("saved",p)
if __name__=="__main__": main()
