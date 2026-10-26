# -*- coding: utf-8 -*-
"""263_slide7_intuitive.py — intuitive story-style page-7 figure (Times New Roman)."""
import os
import numpy as np, pandas as pd, anndata as ad
from scipy.stats import gaussian_kde
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Patch, FancyArrowPatch

OUT=r"D:\桌面\single cell\analysis_output"; STEP1=os.path.join(OUT,"step1")
DEST=os.path.join(OUT,"figures_pub","Slide7_intuitive")
for s in ["svg","tiff","png"]: os.makedirs(os.path.join(DEST,s),exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR",os.path.join(OUT,".mplconfig"))
assert fm.findfont("Times New Roman", fallback_to_default=False)
plt.rcParams.update({"font.family":"serif","font.serif":["Times New Roman"],"font.size":7,
 "axes.titlesize":7.5,"axes.labelsize":6.8,"xtick.labelsize":6.2,"ytick.labelsize":6.2,
 "axes.linewidth":0.6,"svg.fonttype":"none","savefig.facecolor":"white"})

def save(fig,name,dpi=1200):
    fig.savefig(os.path.join(DEST,"svg",name+".svg"),format="svg",bbox_inches="tight")
    fig.savefig(os.path.join(DEST,"tiff",name+".tiff"),format="tiff",dpi=dpi,bbox_inches="tight",pil_kwargs={"compression":"tiff_lzw"})
    fig.savefig(os.path.join(DEST,"png",name+".png"),dpi=220,bbox_inches="tight")
    plt.close(fig); print("saved",name,flush=True)

def letter(ax,s): ax.text(-0.08,1.05,s,transform=ax.transAxes,fontsize=9.5,fontweight="bold",va="top")

# ---------- data ----------
t=ad.read_h5ad(os.path.join(OUT,"adata_hemocyte_traj_geodesic.h5ad"))
ph=pd.read_csv(os.path.join(OUT,"dynamic_subpop_subclusters","PH_subcluster_umap.csv"))
ph.columns=["cell","group","subcluster","umap1","umap2"]
ph["group"]=ph["group"].replace({"Control":"Control","Hemospasia":"Hemospasia"})
meta=ph.set_index("cell")[["group","subcluster"]].rename(columns={"group":"group_ph"})
t.obs=t.obs.join(meta)
phd=t.obs[t.obs.ct=="PH"].copy()
phd["group"]=phd["group_ph"]
phd["subcluster"]=pd.to_numeric(phd["subcluster"],errors="coerce")

# ---------- composite layout ----------
fig=plt.figure(figsize=(7.2,8.4))
gs=fig.add_gridspec(3,3,height_ratios=[1.15,0.95,1.05],hspace=0.55,wspace=0.42,
                    left=0.07,right=0.97,top=0.94,bottom=0.06)
fig.text(0.5,0.975,"PH is the proliferative start; hemospasia shifts PH toward late states",
         ha="center",va="center",fontsize=8.5,fontweight="bold")
fig.text(0.5,0.655,"Hemospasia remodels functional modules and shared DEG programs",
         ha="center",va="center",fontsize=8.5,fontweight="bold")
fig.text(0.5,0.335,"Niche support and cross-organ programme correspondence",
         ha="center",va="center",fontsize=8.5,fontweight="bold")

# a density 2x2
subgs=gs[0,0].subgridspec(2,2,hspace=0.55,wspace=0.35)
for k in range(4):
    ax=fig.add_subplot(subgs[k//2,k%2])
    for grp,col in [("Control","#4C78A8"),("Hemospasia","#E45756")]:
        v=phd.loc[(phd.subcluster==k)&(phd.group==grp),"dpt_pseudotime"].dropna().values
        if len(v)>10:
            x=np.linspace(max(0,v.min()-0.01),min(1,v.max()+0.01),120)
            ax.plot(x,gaussian_kde(v)(x),color=col,lw=1.0)
            ax.axvline(np.median(v),color=col,lw=0.5,ls=":")
    ax.set_title(f"PH_s{k}",fontsize=6.5); ax.set_yticks([])
    ax.set_xlabel("DPT",labelpad=1); ax.tick_params(length=1.5)
    for s in ["top","right","left"]: ax.spines[s].set_visible(False)
    if k==0:
        ax.legend(handles=[Patch(color="#4C78A8",label="Control"),Patch(color="#E45756",label="Hemospasia")],
                  frameon=False,fontsize=5,loc="upper right")
    letter(ax,"a" if k==0 else "")

# b heatmap
h=pd.read_csv(os.path.join(OUT,"Fig5_PH_dpt_heatmap_data.csv"))
piv=h.pivot_table(index="gene",columns="bin",values="zscore")
gene_order=h.groupby("gene")["median_dpt"].mean().sort_values().index.tolist()
piv=piv.reindex(gene_order)
ax=fig.add_subplot(gs[0,1])
im=ax.imshow(piv.values,aspect="auto",cmap="RdBu_r",vmin=-2,vmax=2)
ax.set_yticks(range(len(piv.index))); ax.set_yticklabels(piv.index,fontsize=5.0)
ax.set_xticks([0,len(piv.columns)-1]); ax.set_xticklabels(["early","late"])
ax.set_xlabel("Pseudotime",labelpad=1); ax.set_ylabel("Gene",labelpad=1)
ax.annotate("",xy=(len(piv.columns)-0.5,-0.9),xytext=(-0.5,-0.9),
            arrowprops=dict(arrowstyle="->",lw=0.9,color="#333"),annotation_clip=False)
cb=plt.colorbar(im,ax=ax,fraction=0.035,pad=0.02); cb.set_label("z",fontsize=5.5); cb.ax.tick_params(labelsize=5)
letter(ax,"b")

# c cell cycle
cc=pd.read_csv(os.path.join(OUT,"cellcycle_phase_composition.csv")).set_index("ct").reindex(["PH","GH","MH1","MH2"])
ax=fig.add_subplot(gs[0,2])
left=np.zeros(len(cc))
for phas,col in zip(["G1","S","G2M"],["#BAB0AC","#72B7B2","#E45756"]):
    ax.barh(cc.index,cc[phas],left=left,color=col,label=phas,height=0.55); left+=cc[phas].values
ax.set_xlim(0,100); ax.set_xlabel("Cells (%)",labelpad=1); ax.invert_yaxis()
ax.legend(frameon=False,fontsize=5,ncol=3,loc="lower center",bbox_to_anchor=(0.5,-0.40))
for s in ["top","right"]: ax.spines[s].set_visible(False)
letter(ax,"c")

# d module heatmap
me=pd.read_csv(os.path.join(OUT,"Fig5_coexpression_module_effects.csv"))
piv=me.pivot_table(index="parent",columns="module",values="cohens_d")
ax=fig.add_subplot(gs[1,0:2])
im=ax.imshow(piv.values,aspect="auto",cmap="RdBu_r",vmin=-1.5,vmax=1.5)
ax.set_xticks(range(len(piv.columns))); ax.set_xticklabels(piv.columns)
ax.set_yticks(range(len(piv.index))); ax.set_yticklabels(piv.index)
for i in range(piv.shape[0]):
    for j in range(piv.shape[1]):
        v=piv.values[i,j]
        if np.isfinite(v): ax.text(j,i,f"{v:.1f}",ha="center",va="center",fontsize=5,color="white" if abs(v)>0.9 else "#333")
cb=plt.colorbar(im,ax=ax,fraction=0.025,pad=0.015); cb.set_label("Cohen's d",fontsize=5.5); cb.ax.tick_params(labelsize=5)
letter(ax,"d")

# e DEG sharing
deg=pd.read_csv(os.path.join(OUT,"Fig5_DEG_shared_unique.csv"))
cnt=deg.groupby("combination").size().sort_values(ascending=True).tail(8)
ax=fig.add_subplot(gs[1,2])
ax.barh(cnt.index,cnt.values,color="#4C78A8",height=0.6)
for i,v in enumerate(cnt.values): ax.text(v+1,i,str(v),va="center",fontsize=5.5)
ax.set_xlim(0,cnt.max()*1.22); ax.set_xlabel("DEGs",labelpad=1)
ax.tick_params(axis="y",labelsize=5.0)
for s in ["top","right"]: ax.spines[s].set_visible(False)
letter(ax,"e")

# f bipartite LR
lr=pd.read_csv(os.path.join(STEP1,"step1_lr_communication_sonicparanoid.csv"))
agg=lr.groupby(["source","target"],as_index=False)["prob"].sum()
srcs=["TS","MH2","MyoEpi","M"]; tgts=["PH","GH"]
s_y={x:y for x,y in zip(srcs,np.linspace(1,0,len(srcs)))}; t_y={x:y for x,y in zip(tgts,np.linspace(0.75,0.25,len(tgts)))}
cols={"TS":"#4C78A8","MH2":"#F58518","MyoEpi":"#54A24B","M":"#B279A2"}
ax=fig.add_subplot(gs[2,0]); ax.set_xlim(0,1); ax.set_ylim(-0.1,1.1); ax.axis("off")
for s in srcs:
    ax.scatter(0.12,s_y[s],s=120,color=cols[s],zorder=3)
    ax.text(0.08,s_y[s],s,ha="right",va="center",fontsize=6)
for t_ in tgts:
    ax.scatter(0.88,t_y[t_],s=140,color="#E45756" if t_=="PH" else "#72B7B2",zorder=3)
    ax.text(0.92,t_y[t_],t_,ha="left",va="center",fontsize=6)
for r in agg.itertuples():
    if r.prob<0.01: continue
    p0=np.array([0.14,s_y[r.source]]); p1=np.array([0.86,t_y[r.target]])
    ctrl=np.array([0.5,(p0[1]+p1[1])/2+0.06])
    path=Path([p0,ctrl,p1],[Path.MOVETO,Path.CURVE3,Path.CURVE3])
    ax.add_patch(PathPatch(path,fill=False,color=cols[r.source],lw=0.4+r.prob*7,alpha=0.65))
ax.text(0.5,1.06,"Niche cells",ha="center",fontsize=6.5,fontweight="bold")
ax.text(0.5,0.98,"→ Hemocytes",ha="center",fontsize=6.5)
letter(ax,"f")

# g alluvial
au=pd.read_csv(os.path.join(OUT,"Fig_LyO_mapping_AUROC.csv"))
au=au[(au.auroc>0.70)&(au.support_genes>=2)].copy()
keep=["Progenitor","Macrophage-like 2","Hemocytes","ECM-producing","Lymphocyte precursor-like"]
au=au[au.lyo_name.isin(keep)]
left=sorted(au.lyo_name.unique(),key=lambda x: -au[au.lyo_name==x].auroc.max())
right=[c for c in ["PH","GH","TS","MH2","R","B","MyoEpi","PM-like"] if c in au.target_group.unique()]
ly={n:y for n,y in zip(left,np.linspace(1,0,len(left)))}; ry={n:y for n,y in zip(right,np.linspace(1,0,len(right)))}
lcol=dict(zip(left,["#C44E52","#F58518","#4C78A8","#54A24B","#B279A2"]))
ax=fig.add_subplot(gs[2,1:3]); ax.set_xlim(0,1); ax.set_ylim(-0.05,1.08); ax.axis("off")
for n in left:
    ax.scatter(0.08,ly[n],s=90,color=lcol[n],zorder=3); ax.text(0.045,ly[n],n,ha="right",va="center",fontsize=5.2)
for n in right:
    ax.scatter(0.92,ry[n],s=80,color="#777777",zorder=3); ax.text(0.955,ry[n],n,ha="left",va="center",fontsize=5.2)
for r in au.itertuples():
    p0=np.array([0.10,ly[r.lyo_name]]); p1=np.array([0.90,ry[r.target_group]])
    ctrl=np.array([0.5,(p0[1]+p1[1])/2])
    path=Path([p0,ctrl,p1],[Path.MOVETO,Path.CURVE3,Path.CURVE3])
    ax.add_patch(PathPatch(path,fill=False,color=lcol[r.lyo_name],lw=0.5+(r.auroc-0.5)*4,alpha=0.55))
ax.text(0.08,1.05,"LyO programme",ha="center",fontsize=6.2,fontweight="bold")
ax.text(0.92,1.05,"Hs state",ha="center",fontsize=6.2,fontweight="bold")
letter(ax,"g")

fig.savefig(os.path.join(DEST,"png","Fig7_intuitive_composite.png"),dpi=300,bbox_inches="tight")
fig.savefig(os.path.join(DEST,"svg","Fig7_intuitive_composite.svg"),format="svg",bbox_inches="tight")
fig.savefig(os.path.join(DEST,"tiff","Fig7_intuitive_composite.tiff"),format="tiff",dpi=1200,bbox_inches="tight",pil_kwargs={"compression":"tiff_lzw"})
print("composite saved")
