# -*- coding: utf-8 -*-
"""201_bulk_deconvolution_full.py
Deconvolve the complete GHRL25030101-2 bulk matrix against the 10 single-cell groups.
"""
import os, sys, re
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
from scipy.optimize import nnls
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
OUT=r"D:\桌面\single cell\analysis_output"
PREV=os.path.join(OUT,"figures_pub","sc_summary_new_names","individual_previews")
os.makedirs(PREV,exist_ok=True)
GROUPS={"0":"TS","1":"GH","2":"PH","3":"MH2","4":"PM-like","5":"F","6":"R","7":"B","8":"MyoEpi","9":"M"}
ORDER=["PH","GH","TS","MH2","F","R","B","M","MyoEpi","PM-like"]

def main():
    B=pd.read_csv(os.path.join(OUT,"bulk_GHRL25030101","5.ExpressionStat","all.genes.expression.annot.xls"),sep="\t",usecols=["Symbol","haemospasia_count","untreated_count"])
    B["Symbol"]=B.Symbol.astype(str); B=B[B.Symbol!="--"].copy()
    B=B.groupby("Symbol",as_index=False).agg(untreated_count=("untreated_count","sum"),haemospasia_count=("haemospasia_count","sum"))
    A=ad.read_h5ad(os.path.join(OUT,"adata_soupx.h5ad")); A.X=A.layers["counts"].copy(); sc.pp.normalize_total(A,target_sum=1e4); sc.pp.log1p(A)
    ct=A.obs.vendor_cluster.astype(str).map(GROUPS).values
    sym2gid={}
    for gid,s in zip(A.var_names.astype(str),A.var["sym"].astype(str)): sym2gid.setdefault(s,gid)
    # aggregate single-cell by symbol
    syms=sorted(set(B.Symbol) & set(sym2gid))
    print("bulk symbols",B.Symbol.nunique(),"single symbols",len(sym2gid),"matched",len(syms))
    X=np.asarray(A[:,[sym2gid[s] for s in syms]].X.toarray(),dtype=np.float32)
    S=[]
    for g in ORDER: S.append(X[ct==g].mean(0))
    S=np.vstack(S).T
    S=S/(S.sum(0,keepdims=True)+1e-9)
    b=B.set_index("Symbol").loc[syms]
    b0=b.untreated_count.values.astype(float); b6=b.haemospasia_count.values.astype(float)
    b0=b0/(b0.sum()+1e-9); b6=b6/(b6.sum()+1e-9)
    def fit(b):
        coef,_=nnls(S,b); return coef/(coef.sum()+1e-9)
    p0=fit(b0); p6=fit(b6)
    # marker score: top 10 markers/symbol
    M=pd.read_csv(os.path.join(OUT,"top_markers_10groups.csv")); M=M[M["rank"]<=10]
    rows=[]
    for g in ORDER:
        mg=sorted(set(M[M.grp==g].sym.astype(str)) & set(B.Symbol))
        d=b.loc[mg] if mg else b.iloc[0:0]
        if len(mg)>=3:
            s0=float(np.log1p(d.untreated_count).mean()); s6=float(np.log1p(d.haemospasia_count).mean())
            n0=float(np.log1p(d.untreated_count).mean()/(np.log1p(b.untreated_count).mean()+1e-9))
            n6=float(np.log1p(d.haemospasia_count).mean()/(np.log1p(b.haemospasia_count).mean()+1e-9))
        else: s0=s6=n0=n6=np.nan
        rows.append({"group":g,"nnls_untreated":p0[ORDER.index(g)],"nnls_haemospasia":p6[ORDER.index(g)],"n_markers":len(mg),"marker_score_untreated":s0,"marker_score_haemospasia":s6})
    R=pd.DataFrame(rows); R.to_csv(os.path.join(OUT,"bulk_deconvolution_full_results.csv"),index=False,encoding="utf-8-sig")
    print(R.round(4).to_string(index=False))
    fig=plt.figure(figsize=(5.3,3.4),facecolor="white")
    ax=fig.add_axes([0.20,0.20,0.77,0.70]); y=np.arange(len(R))
    ax.barh(y-0.18,R.nnls_untreated,height=0.34,color="#4C78A8",label="Untreated")
    ax.barh(y+0.18,R.nnls_haemospasia,height=0.34,color="#D6604D",label="Haemospasia")
    ax.set_yticks(y); ax.set_yticklabels(R.group,fontsize=5.6); ax.invert_yaxis()
    ax.set_xlabel("Estimated proportion (NNLS)",fontsize=5.8)
    ax.legend(frameon=False,fontsize=5.0,loc="lower right")
    ax.set_title("Bulk deconvolution: untreated vs haemospasia",fontsize=6.8,pad=4)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)
    ax.tick_params(axis="x",labelsize=4.8,length=1.0); ax.tick_params(axis="y",length=0)
    p=os.path.join(PREV,"06b_bulk_deconvolution_full_preview.png"); fig.savefig(p,dpi=400,bbox_inches="tight",pad_inches=0.04,facecolor="white"); plt.close(fig); print("saved",p)

if __name__=="__main__": main()
