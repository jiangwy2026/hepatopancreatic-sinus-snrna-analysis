# -*- coding: utf-8 -*-
"""198_prior_pathway_modules_preview.py
Map prior protein/metabolite programmes and betaine pathway to single-cell groups.
"""
import os, sys
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
OUT=r"D:\桌面\single cell\analysis_output"; PREV=os.path.join(OUT,"figures_pub","sc_summary_new_names","individual_previews")
GROUPS={"0":"TS","1":"GH","2":"PH","3":"MH2","4":"PM-like","5":"F","6":"R","7":"B","8":"MyoEpi","9":"M"}
ORDER=["PH","GH","TS","MH2","F","R","B","M","MyoEpi","PM-like"]
MODULES={
 "Protein synthesis\n& cycle initiation":["ncbi_113812494","ncbi_113811156","ncbi_113819376","ncbi_113823483","ncbi_113820946","ncbi_113821402","ncbi_113822791","ncbi_113804001"],
 "Adhesion\ncytoskeleton":["ncbi_113806735","ncbi_113813659"],
 "Antioxidant\n/ iron":["ncbi_113823550","ncbi_113830370","ncbi_113821440","ncbi_113823144"],
 "Betaine / choline":["ncbi_113819630","ncbi_113817275","ncbi_113805144","ncbi_113808193","ncbi_113806314","ncbi_113826080"],
}
BETAINE={"bhmt":["ncbi_113819630"],"Aldh9a1":["ncbi_113817275"],"Bbox1":["ncbi_113805144"],"ChT":["ncbi_113808193"],"Nfat5":["ncbi_113806314"],"CHDH":["ncbi_113826080"]}

def main():
    A=ad.read_h5ad(os.path.join(OUT,"adata_soupx.h5ad")); A.X=A.layers["counts"].copy(); sc.pp.normalize_total(A,target_sum=1e4); sc.pp.log1p(A)
    ct=A.obs.vendor_cluster.astype(str).map(GROUPS).values; grp=A.obs["sample"].astype(str).values
    gid2i={str(g):i for i,g in enumerate(A.var_names.astype(str))}
    X=np.asarray(A.X.toarray(),dtype=np.float32)
    # z-score genes then average within module
    scores={}
    for name,ids in MODULES.items():
        ids=[g for g in ids if g in gid2i]
        Z=(X[:,[gid2i[g] for g in ids]]-X[:,[gid2i[g] for g in ids]].mean(0))/(X[:,[gid2i[g] for g in ids]].std(0)+1e-6)
        scores[name]=Z.mean(1)
    S=pd.DataFrame(scores); S["group"]=ct; S["sample"]=grp
    M=S.groupby("group")[list(MODULES)].mean().reindex(ORDER)
    M.to_csv(os.path.join(OUT,"prior_program_module_scores_by_group.csv"),encoding="utf-8-sig")
    fig=plt.figure(figsize=(7.2,4.5),facecolor="white")
    gs=GridSpec(2,2,height_ratios=[1.05,1.0],width_ratios=[1.0,1.0],left=0.10,right=0.97,top=0.90,bottom=0.10,wspace=0.42,hspace=0.55)
    axa=fig.add_subplot(gs[0,0])
    Z=(M-M.mean(0))/ (M.std(0)+1e-9)
    im=axa.imshow(Z.T.values,aspect="auto",cmap="RdBu_r",vmin=-2,vmax=2)
    axa.set_yticks(range(len(MODULES))); axa.set_yticklabels(list(MODULES),fontsize=4.8)
    axa.set_xticks(range(len(ORDER))); axa.set_xticklabels(ORDER,fontsize=5.0,rotation=35,ha="right")
    axa.set_title("Prior protein/metabolite programmes",fontsize=6.2,pad=3)
    axa.tick_params(length=0)
    for sp in axa.spines.values(): sp.set_visible(False)
    cb=fig.colorbar(im,ax=axa,fraction=0.035,pad=0.02); cb.set_label("row z-score",fontsize=4.6); cb.ax.tick_params(labelsize=4.3,length=0.6); cb.outline.set_visible(False)
    # PH module score by sample
    axb=fig.add_subplot(gs[0,1])
    ph=S[S.group=="PH"].groupby("sample")[list(MODULES)].mean().reindex(["blank","treated"])
    xx=np.arange(len(MODULES)); axb.bar(xx-0.18,ph.loc["blank"],width=0.34,color="#4C78A8",label="Control")
    axb.bar(xx+0.18,ph.loc["treated"],width=0.34,color="#D6604D",label="Hemospasia")
    axb.set_xticks(xx); axb.set_xticklabels(["Protein\nsynthesis","Adhesion","Antioxidant","Betaine"],fontsize=4.8)
    axb.set_ylabel("Module score",fontsize=4.8); axb.tick_params(axis="y",labelsize=4.3,length=0.8)
    axb.legend(frameon=False,fontsize=4.5,loc="upper left"); axb.set_title("PH response",fontsize=6.0,pad=3)
    for sp in ("top","right"): axb.spines[sp].set_visible(False)
    # betaine dot plot
    axc=fig.add_subplot(gs[1,:])
    rows=[]; genes=list(BETAINE)
    for g in genes:
        gid=BETAINE[g][0]
        v=X[:,gid2i[gid]]
        for c in ORDER:
            m=ct==c
            rows.append({"gene":g,"group":c,"mean":float(v[m].mean()),"pct":float((v[m]>0).mean()*100)})
    B=pd.DataFrame(rows)
    B.to_csv(os.path.join(OUT,"betaine_pathway_expression.csv"),index=False,encoding="utf-8-sig")
    means=B.pivot(index="gene",columns="group",values="mean").reindex(genes)[ORDER]
    pcts=B.pivot(index="gene",columns="group",values="pct").reindex(genes)[ORDER]
    Z=(means-means.mean(1).values[:,None])/(means.std(1).values[:,None]+1e-9)
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            axc.scatter(j,i,s=max(2,pcts.iloc[i,j]*2.0),c=Z.iloc[i,j],cmap="RdBu_r",vmin=-2,vmax=2,edgecolors="white",linewidths=0.3)
    axc.set_xticks(range(len(ORDER))); axc.set_xticklabels(ORDER,fontsize=5.0,rotation=35,ha="right")
    axc.set_yticks(range(len(genes))); axc.set_yticklabels(genes,fontsize=5.0)
    axc.set_xlim(-0.5,len(ORDER)-0.5); axc.set_ylim(len(genes)-0.5,-0.5); axc.tick_params(length=0)
    for sp in axc.spines.values(): sp.set_visible(False)
    axc.set_title("Betaine / choline pathway",fontsize=6.2,pad=4)
    fig.text(0.50,0.025,"Bubble size = % detected; color = row-scaled mean expression.",ha="center",va="center",fontsize=4.8,color="#555555")
    p=os.path.join(PREV,"09_prior_pathways_betaine_preview.png"); fig.savefig(p,dpi=400,bbox_inches="tight",pad_inches=0.04,facecolor="white"); plt.close(fig); print("saved",p)

if __name__=="__main__": main()
