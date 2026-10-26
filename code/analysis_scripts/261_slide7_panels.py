# -*- coding: utf-8 -*-
"""261_slide7_panels.py — re-analysis panels for PPT page 7 (Times New Roman, SVG + 1200 dpi TIFF)."""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import Patch

OUT=r"D:\桌面\single cell\analysis_output"; STEP1=os.path.join(OUT,"step1")
DEST=os.path.join(OUT,"figures_pub","Slide7_panels")
for s in ["svg","tiff","png"]: os.makedirs(os.path.join(DEST,s),exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR",os.path.join(OUT,".mplconfig"))
assert fm.findfont("Times New Roman", fallback_to_default=False)
plt.rcParams.update({"font.family":"serif","font.serif":["Times New Roman"],"font.size":7,
 "axes.titlesize":8,"axes.labelsize":7,"xtick.labelsize":6.5,"ytick.labelsize":6.5,
 "axes.linewidth":0.6,"svg.fonttype":"none","savefig.facecolor":"white"})

def save(fig,name):
    fig.savefig(os.path.join(DEST,"svg",name+".svg"),format="svg",bbox_inches="tight")
    fig.savefig(os.path.join(DEST,"tiff",name+".tiff"),format="tiff",dpi=1200,bbox_inches="tight",pil_kwargs={"compression":"tiff_lzw"})
    fig.savefig(os.path.join(DEST,"png",name+".png"),dpi=200,bbox_inches="tight")
    plt.close(fig); print("saved",name,flush=True)

def letter(ax,s): ax.text(-0.10,1.04,s,transform=ax.transAxes,fontsize=10,fontweight="bold",va="top")

# a PH pseudotime by subcluster
ps=pd.read_csv(os.path.join(OUT,"Fig5_dynamic_subcluster_pseudotime_summary.csv"))
ps=ps[ps.parent=="PH"].copy()
order=["PH_s1","PH_s0","PH_s2","PH_s3"]
fig,ax=plt.subplots(figsize=(3.8,2.8))
for i,k in enumerate(order):
    for j,(grp,col) in enumerate([("Control","#4C78A8"),("Hemospasia","#E45756")]):
        r=ps[(ps.key==k)&(ps.group==grp)]
        if len(r)==0: continue
        r=r.iloc[0]; y=i+(j-0.5)*0.22
        ax.errorbar(r.median_pseudotime,y,xerr=[[r.median_pseudotime-r.q25],[r.q75-r.median_pseudotime]],
                    fmt="o",ms=3.2,color=col,capsize=1.5,lw=0.8)
ax.set_yticks(range(len(order))); ax.set_yticklabels(order); ax.invert_yaxis()
ax.set_xlabel("DPT pseudotime"); ax.set_ylabel("PH substate")
ax.legend(handles=[Patch(color="#4C78A8",label="Control"),Patch(color="#E45756",label="Hemospasia")],
          frameon=False,loc="lower right",fontsize=6)
for s in ["top","right"]: ax.spines[s].set_visible(False)
letter(ax,"a"); fig.tight_layout(); save(fig,"panel_a_PH_pseudotime")

# b PH pseudotime heatmap
h=pd.read_csv(os.path.join(OUT,"Fig5_PH_dpt_heatmap_data.csv"))
piv=h.pivot_table(index="gene",columns="bin",values="zscore")
gene_order=h.groupby("gene")["median_dpt"].mean().sort_values().index.tolist()
piv=piv.reindex(gene_order)
fig,ax=plt.subplots(figsize=(3.8,3.2))
im=ax.imshow(piv.values,aspect="auto",cmap="RdBu_r",vmin=-2,vmax=2)
ax.set_yticks(range(len(piv.index))); ax.set_yticklabels(piv.index,fontsize=5.8)
ax.set_xticks([0,len(piv.columns)-1]); ax.set_xticklabels(["early","late"])
ax.set_xlabel("Pseudotime"); ax.set_ylabel("Gene")
cb=plt.colorbar(im,ax=ax,fraction=0.04,pad=0.02); cb.set_label("z-score",fontsize=6); cb.ax.tick_params(labelsize=5.5)
letter(ax,"b"); fig.tight_layout(); save(fig,"panel_b_PH_heatmap")

# c cell cycle
cc=pd.read_csv(os.path.join(OUT,"cellcycle_phase_composition.csv"))
cc=cc.set_index("ct").reindex(["PH","GH","MH1","MH2"])
fig,ax=plt.subplots(figsize=(3.2,2.6))
left=np.zeros(len(cc))
for ph,col in zip(["G1","S","G2M"],["#BAB0AC","#72B7B2","#E45756"]):
    ax.barh(cc.index,cc[ph],left=left,color=col,label=ph,height=0.6); left+=cc[ph].values
ax.set_xlabel("Cells per phase (%)"); ax.set_xlim(0,100)
ax.legend(frameon=False,fontsize=6,ncol=3,loc="lower center",bbox_to_anchor=(0.5,-0.42))
for s in ["top","right"]: ax.spines[s].set_visible(False)
letter(ax,"c"); fig.tight_layout(); save(fig,"panel_c_cellcycle")

# d module effects
me=pd.read_csv(os.path.join(OUT,"Fig5_coexpression_module_effects.csv"))
piv=me.pivot_table(index="parent",columns="module",values="cohens_d")
fig,ax=plt.subplots(figsize=(3.6,2.6))
im=ax.imshow(piv.values,aspect="auto",cmap="RdBu_r",vmin=-1.5,vmax=1.5)
ax.set_xticks(range(len(piv.columns))); ax.set_xticklabels(piv.columns)
ax.set_yticks(range(len(piv.index))); ax.set_yticklabels(piv.index)
for i in range(piv.shape[0]):
    for j in range(piv.shape[1]):
        v=piv.values[i,j]
        if np.isfinite(v): ax.text(j,i,f"{v:.1f}",ha="center",va="center",fontsize=5.5,color="white" if abs(v)>0.9 else "#333")
cb=plt.colorbar(im,ax=ax,fraction=0.04,pad=0.02); cb.set_label("Cohen's d",fontsize=6); cb.ax.tick_params(labelsize=5.5)
letter(ax,"d"); fig.tight_layout(); save(fig,"panel_d_modules")

# e DEG sharing
deg=pd.read_csv(os.path.join(OUT,"Fig5_DEG_shared_unique.csv"))
cnt=deg.groupby("combination").size().sort_values(ascending=True)
fig,ax=plt.subplots(figsize=(3.6,2.8))
ax.barh(cnt.index,cnt.values,color="#4C78A8",height=0.65)
for i,v in enumerate(cnt.values): ax.text(v+1,i,str(v),va="center",fontsize=6)
ax.set_xlabel("Number of DEGs"); ax.set_xlim(0,cnt.max()*1.18)
for s in ["top","right"]: ax.spines[s].set_visible(False)
letter(ax,"e"); fig.tight_layout(); save(fig,"panel_e_DEG_sharing")

# f validated LR
lr=pd.read_csv(os.path.join(STEP1,"step1_lr_communication_sonicparanoid.csv"))
sig=lr[lr.p_perm<0.05].copy()
sig["label"]=sig["pair"]+"  ("+sig["source"]+"→"+sig["target"]+")"
top=sig.sort_values("prob",ascending=False).drop_duplicates("label").head(10).sort_values("prob")
fig,ax=plt.subplots(figsize=(4.0,3.0))
ax.barh(top["label"],top["prob"],color="#C44E52",height=0.65)
for i,(v,p) in enumerate(zip(top["prob"],top["p_perm"])):
    ax.text(v+0.003,i,f"P={p:.1e}",va="center",fontsize=5.5)
ax.set_xlabel("Communication probability"); ax.set_xlim(0,top["prob"].max()*1.35)
for s in ["top","right"]: ax.spines[s].set_visible(False)
letter(ax,"f"); fig.tight_layout(); save(fig,"panel_f_validated_LR")

# g LyO AUROC heatmap
au=pd.read_csv(os.path.join(OUT,"Fig_LyO_mapping_AUROC.csv"))
keep=["ECM-producing","Macrophage-like 1","Macrophage-like 2","Lymphocyte precursor-like","Capsular","Hemocytes","SLC-rich","Epithelial","Progenitor"]
au=au[au.lyo_name.isin(keep)]
piv=au.pivot_table(index="lyo_name",columns="target_group",values="auroc").reindex(keep)
tg=[c for c in ["PH","GH","TS","MH2","F","R","B","M","MyoEpi","PM-like"] if c in piv.columns]
piv=piv[tg]
fig,ax=plt.subplots(figsize=(3.8,3.0))
im=ax.imshow(piv.values,aspect="auto",cmap="viridis",vmin=0.5,vmax=1.0)
ax.set_xticks(range(len(tg))); ax.set_xticklabels(tg,rotation=60)
ax.set_yticks(range(len(piv.index))); ax.set_yticklabels(piv.index,fontsize=5.8)
cb=plt.colorbar(im,ax=ax,fraction=0.04,pad=0.02); cb.set_label("AUROC",fontsize=6); cb.ax.tick_params(labelsize=5.5)
letter(ax,"g"); fig.tight_layout(); save(fig,"panel_g_LyO_mapping")

print("Slide 7 panels done ->",DEST)
