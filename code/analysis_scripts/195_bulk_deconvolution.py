# -*- coding: utf-8 -*-
"""195_bulk_deconvolution.py
Use the prior Hs bulk TPM table (Supplemental Table 5, tpm>100)
to estimate 10-cluster proportions in untreated vs haemospasia.
"""
import os, sys
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
from scipy.optimize import nnls
from scipy.stats import zscore
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)

ROOT=r"D:\桌面\single cell"; OUT=os.path.join(ROOT,"analysis_output")
PREV=os.path.join(OUT,"figures_pub","sc_summary_new_names","individual_previews")
os.makedirs(PREV,exist_ok=True)
GROUPS={"0":"TS","1":"GH","2":"PH","3":"MH2","4":"PM-like","5":"F","6":"R","7":"B","8":"MyoEpi","9":"M"}
ORDER=["PH","GH","TS","MH2","F","R","B","M","MyoEpi","PM-like"]
COL={"PH":"#B2182B","GH":"#E08214","TS":"#D6604D","MH2":"#8C510A","F":"#006D2C","R":"#41AB5D","B":"#78C679","M":"#2166AC","MyoEpi":"#92C5DE","PM-like":"#7B3294"}

def main():
    bulk=pd.read_csv(os.path.join(OUT,"prior_bulk_tpm_genes.csv"))
    bulk["gene_id"]=bulk.gene_id.astype(str)
    A=ad.read_h5ad(os.path.join(OUT,"adata_soupx.h5ad"))
    ct=A.obs.vendor_cluster.astype(str).map(GROUPS).values
    A.X=A.layers["counts"].copy(); sc.pp.normalize_total(A,target_sum=1e4); sc.pp.log1p(A)
    gid2i={str(g):i for i,g in enumerate(A.var_names.astype(str))}
    genes=[g for g in bulk.gene_id if g in gid2i]
    bulk=bulk[bulk.gene_id.isin(genes)].drop_duplicates("gene_id").set_index("gene_id")
    X=np.asarray(A[:,genes].X.toarray(),dtype=float)
    S=[]
    for g in ORDER:
        S.append(X[ct==g].mean(0))
    S=np.vstack(S).T  # genes x clusters
    # normalize columns to relative gene expression
    S=S/(S.sum(0,keepdims=True)+1e-9)
    b0=bulk.loc[genes,"untreated_tpm"].values.astype(float); b6=bulk.loc[genes,"haemospasia_tpm"].values.astype(float)
    b0=b0/(b0.sum()+1e-9); b6=b6/(b6.sum()+1e-9)
    def fit(b):
        coef,_=nnls(S,b)
        if coef.sum()<=0: return coef
        return coef/coef.sum()
    p0=fit(b0); p6=fit(b6)
    # marker score: mean z-scored bulk TPM across top markers
    M=pd.read_csv(os.path.join(OUT,"top_markers_10groups.csv"))
    M=M[M["rank"]<=10].copy()
    sym2gid={}
    for gid,s in zip(A.var_names.astype(str),A.var["sym"].astype(str)): sym2gid.setdefault(s,gid)
    # use gid from marker table directly
    rows=[]
    for g in ORDER:
        mg=M[M.grp==g]["gid"].astype(str).tolist()
        mg=[x for x in mg if x in bulk.index]
        if len(mg)>=3:
            v=bulk.loc[mg,["untreated_tpm","haemospasia_tpm"]].mean(1)
            # simple relative marker score
            ms0=float(np.log1p(bulk.loc[mg,"untreated_tpm"]).mean())
            ms6=float(np.log1p(bulk.loc[mg,"haemospasia_tpm"]).mean())
        else:
            ms0=ms6=np.nan
        rows.append({"group":g,"nnls_untreated":p0[ORDER.index(g)],"nnls_haemospasia":p6[ORDER.index(g)],
                     "marker_logmean_untreated":ms0,"marker_logmean_haemospasia":ms6})
    R=pd.DataFrame(rows)
    R.to_csv(os.path.join(OUT,"bulk_deconvolution_results.csv"),index=False,encoding="utf-8-sig")
    print(R.round(4).to_string(index=False))
    print("genes used",len(genes))
    fig=plt.figure(figsize=(5.2,3.4),facecolor="white")
    ax=fig.add_axes([0.20,0.20,0.77,0.70])
    y=np.arange(len(R))
    ax.barh(y-0.18,R.nnls_untreated,height=0.34,color="#4C78A8",label="Untreated")
    ax.barh(y+0.18,R.nnls_haemospasia,height=0.34,color="#D6604D",label="Haemospasia")
    ax.set_yticks(y); ax.set_yticklabels(R.group,fontsize=5.8); ax.invert_yaxis()
    ax.set_xlabel("Estimated proportion (NNLS)",fontsize=5.8)
    ax.legend(frameon=False,fontsize=5.2,loc="lower right")
    ax.set_title("Prior Hs bulk deconvolution",fontsize=7.0,pad=5)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)
    ax.tick_params(axis="x",labelsize=5.0,length=1.1); ax.tick_params(axis="y",length=0)
    p=os.path.join(PREV,"06_bulk_deconvolution_preview.png")
    fig.savefig(p,dpi=400,bbox_inches="tight",pad_inches=0.04,facecolor="white"); plt.close(fig)
    print("saved",p)

if __name__=="__main__": main()
