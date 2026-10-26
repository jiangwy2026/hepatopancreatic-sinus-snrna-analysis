# -*- coding: utf-8 -*-
"""197_table1_conserved_genes_preview.py - map Table 1 orthologs to the 10 single-cell groups."""
import os, sys, re
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
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
    M=pd.read_csv(os.path.join(OUT,"Table1_ortholog_mapping_full.csv"))
    M=M[M.gene_id.notna()].copy()
    def first_id(x):
        m=re.search(r'(ncbi_\d+)',str(x)); return m.group(1) if m else None
    M["gid"]=M.gene_id.map(first_id)
    M=M[M.gid.notna()].drop_duplicates(["human_gene","human_cluster","shrimp_gene","gid"])
    A=ad.read_h5ad(os.path.join(OUT,"adata_soupx.h5ad")); A.X=A.layers["counts"].copy(); sc.pp.normalize_total(A,target_sum=1e4); sc.pp.log1p(A)
    ct=A.obs.vendor_cluster.astype(str).map(GROUPS).values
    gid2i={str(g):i for i,g in enumerate(A.var_names.astype(str))}
    M=M[M.gid.isin(gid2i)].copy()
    print("mapped rows",len(M),"genes",M.human_gene.nunique())
    rows=[]
    for r in M.itertuples():
        v=A[:,gid2i[r.gid]].X.toarray().ravel()
        for g in ORDER:
            m=ct==g
            rows.append({"human_gene":r.human_gene,"human_cluster":r.human_cluster,"shrimp_gene":r.shrimp_gene,
                         "gid":r.gid,"group":g,"mean":float(v[m].mean()),"pct":float((v[m]>0).mean()*100)})
    R=pd.DataFrame(rows)
    # collapse duplicate human_gene rows by max mean
    R=R.groupby(["human_gene","human_cluster","group"]).agg(mean=("mean","max"),pct=("pct","max")).reset_index()
    R.to_csv(os.path.join(OUT,"Table1_conserved_gene_expression.csv"),index=False,encoding="utf-8-sig")
    genes=R[["human_gene","human_cluster"]].drop_duplicates().sort_values("human_gene").values.tolist()
    fig=plt.figure(figsize=(6.8,4.4),facecolor="white")
    ax=fig.add_axes([0.27,0.16,0.66,0.74])
    means=R.pivot(index=["human_gene","human_cluster"],columns="group",values="mean").reindex([tuple(x) for x in genes])[ORDER]
    pcts=R.pivot(index=["human_gene","human_cluster"],columns="group",values="pct").reindex([tuple(x) for x in genes])[ORDER]
    Z=means.sub(means.mean(axis=1),axis=0).div(means.std(axis=1)+1e-6,axis=0)
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            ax.scatter(j,i,s=max(2.0,pcts.iloc[i,j]*1.8),c=Z.iloc[i,j],cmap="RdBu_r",vmin=-2,vmax=2,edgecolors="white",linewidths=0.3)
    ax.set_xticks(range(len(ORDER))); ax.set_xticklabels(ORDER,fontsize=5.2,rotation=35,ha="right")
    ax.set_yticks(range(len(genes))); ax.set_yticklabels(["%s [%s]"%(g,c) for g,c in genes],fontsize=4.8)
    ax.set_xlim(-0.5,len(ORDER)-0.5); ax.set_ylim(len(genes)-0.5,-0.5)
    ax.set_title("Prior Table 1 conserved genes in single-cell groups",fontsize=6.7,pad=5)
    ax.tick_params(length=0)
    for sp in ax.spines.values(): sp.set_visible(False)
    # colorbar
    cax=fig.add_axes([0.945,0.35,0.014,0.30]); cb=fig.colorbar(plt.cm.ScalarMappable(cmap="RdBu_r",norm=plt.Normalize(-2,2)),cax=cax)
    cb.set_label("row z-score",fontsize=4.8); cb.ax.tick_params(labelsize=4.5,length=0.7); cb.outline.set_visible(False)
    fig.text(0.50,0.035,"Bubble size = % detected; only genes with an ortholog mapping are shown.",ha="center",va="center",fontsize=4.8,color="#555555")
    p=os.path.join(PREV,"08_Table1_conserved_genes_preview.png"); fig.savefig(p,dpi=400,bbox_inches="tight",pad_inches=0.04,facecolor="white"); plt.close(fig); print("saved",p)

if __name__=="__main__": main()
