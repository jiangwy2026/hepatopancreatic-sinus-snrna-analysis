# -*- coding: utf-8 -*-
"""253_make_summary_figure.py
Composite figure supporting the conclusion:
PH shares the lymphoid progenitor program; niche cells do not enter it and
instead provide ECM-integrin support.
"""
import os, sys
import numpy as np, pandas as pd, anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

OUT = r"D:\桌面\single cell\analysis_output"
STEP1 = os.path.join(OUT,"step1"); STEP2 = os.path.join(OUT,"step2"); STEP3 = os.path.join(OUT,"step3")
os.environ.setdefault("MPLCONFIGDIR", os.path.join(OUT,".mplconfig"))
plt.rcParams.update({"font.size":7,"axes.titlesize":8,"axes.labelsize":7,
                     "xtick.labelsize":6.5,"ytick.labelsize":6.5,"pdf.fonttype":42})

a = ad.read_h5ad(os.path.join(STEP3,"step3_lyo_scores.h5ad"))
uc = pd.read_csv(os.path.join(STEP3,"step3_ucell_by_celltype.csv"))
ov = pd.read_csv(os.path.join(STEP3,"step3_hypergeom_overlaps.csv"))
lr = pd.read_csv(os.path.join(STEP1,"step1_lr_communication_sonicparanoid.csv"))
order = ["PH","GH","TS","MH2","MyoEpi","M","F","R","B","PM-like"]
colors = {"PH":"#E45756","GH":"#72B7B2","TS":"#4C78A8","MH2":"#F58518","MyoEpi":"#54A24B",
          "M":"#B279A2","F":"#FF9DA6","R":"#9D755D","B":"#BAB0AC","PM-like":"#79706E"}

fig = plt.figure(figsize=(7.2,7.0))
gs = fig.add_gridspec(3,2, hspace=0.42, wspace=0.28, left=0.08, right=0.97, top=0.95, bottom=0.07)

# A cell type UMAP
ax=fig.add_subplot(gs[0,0])
for ct in order:
    m=(a.obs.cell_type==ct).values
    ax.scatter(a.obsm["X_umap"][m,0],a.obsm["X_umap"][m,1],s=1.2,c=colors[ct],linewidths=0,rasterized=True)
    if m.sum()>20:
        ax.text(np.median(a.obsm["X_umap"][m,0]),np.median(a.obsm["X_umap"][m,1]),ct,fontsize=6,ha="center",
                bbox=dict(boxstyle="round,pad=0.12",fc="white",ec="none",alpha=0.65))
ax.set_title("a  Cell types"); ax.set_xticks([]); ax.set_yticks([])
for s in ax.spines.values(): s.set_visible(False)

# B progenitor UCell UMAP
ax=fig.add_subplot(gs[0,1])
v=a.obs["ucell_LyO_Progenitor"].values
sc=ax.scatter(a.obsm["X_umap"][:,0],a.obsm["X_umap"][:,1],c=v,cmap="viridis",s=1.2,linewidths=0,rasterized=True,
              vmin=np.quantile(v,0.02),vmax=np.quantile(v,0.98))
ax.set_title("b  LyO progenitor program (UCell)"); ax.set_xticks([]); ax.set_yticks([])
for s in ax.spines.values(): s.set_visible(False)
plt.colorbar(sc,ax=ax,fraction=0.042,pad=0.01)

# C boxplot progenitor score
ax=fig.add_subplot(gs[1,0])
data=[a.obs.loc[a.obs.cell_type==c,"ucell_LyO_Progenitor"].values for c in order]
bp=ax.boxplot(data,tick_labels=order,showfliers=False,patch_artist=True,widths=0.62)
for patch,c in zip(bp["boxes"],order): patch.set_facecolor(colors[c]); patch.set_alpha(0.75)
ax.set_title("c  Progenitor program by cell type"); ax.tick_params(axis="x",rotation=60)
ax.set_ylabel("UCell score")

# D hypergeometric fold
ax=fig.add_subplot(gs[1,1])
sub=ov[ov.lyo_set=="LyO_Progenitor"].set_index("compared_to").reindex(["PH_markers","niche_markers","TS_markers","MH2_markers","MyoEpi_markers","M_markers"])
y=np.arange(len(sub)); vals=sub["fold"].fillna(0).values
ax.barh(y,vals,color=["#C44E52" if x>1 else "#BBBBBB" for x in vals])
ax.set_yticks(y); ax.set_yticklabels(["PH","niche","TS","MH2","MyoEpi","M"]); ax.invert_yaxis()
for i,(val,p,k) in enumerate(zip(vals,sub["p_hypergeom"].values,sub["overlap"].values)):
    ax.text(max(val,0.02)+2,i,f"{int(k)} genes, {val:.1f}x\np={p:.1e}",va="center",fontsize=5.5)
ax.set_xlim(0,max(12,vals.max()*1.6)); ax.set_xlabel("Fold enrichment")
ax.set_title("d  Progenitor top50 overlap"); 
for s in ["top","right"]: ax.spines[s].set_visible(False)

# E LR heatmap
ax=fig.add_subplot(gs[2,0])
piv=lr.groupby(["source","target"])["prob"].sum().unstack(fill_value=0).reindex(index=["TS","MH2","MyoEpi","M"],columns=["PH","GH"],fill_value=0)
im=ax.imshow(piv.values,cmap="Reds",aspect="auto")
ax.set_xticks([0,1]); ax.set_xticklabels(["PH","GH"]); ax.set_yticks(range(4)); ax.set_yticklabels(["TS","MH2","MyoEpi","M"])
for i in range(4):
    for j in range(2):
        ax.text(j,i,f"{piv.values[i,j]:.2f}",ha="center",va="center",fontsize=6,
                color="white" if piv.values[i,j]>piv.values.max()*0.6 else "#333333")
ax.set_title("e  Niche -> hemocyte ECM-integrin"); plt.colorbar(im,ax=ax,fraction=0.042,pad=0.01)

# F projected pseudotime
ax=fig.add_subplot(gs[2,1])
data=[a.obs.loc[a.obs.cell_type==c,"projected_pseudotime"].values for c in order]
bp=ax.boxplot(data,tick_labels=order,showfliers=False,patch_artist=True,widths=0.62)
for patch,c in zip(bp["boxes"],order): patch.set_facecolor(colors[c]); patch.set_alpha(0.75)
ax.set_title("f  Projection onto PH->GH trajectory"); ax.tick_params(axis="x",rotation=60)
ax.set_ylabel("Projected pseudotime")
ax.annotate("PH start",xy=(1,0.27),xytext=(1.5,0.22),fontsize=5.5,arrowprops=dict(arrowstyle="->",lw=0.6))
ax.annotate("GH end",xy=(2,0.40),xytext=(2.6,0.45),fontsize=5.5,arrowprops=dict(arrowstyle="->",lw=0.6))

fig.savefig(os.path.join(OUT,"figures_pub","Fig_niche_model_summary.png"),dpi=600,bbox_inches="tight")
fig.savefig(os.path.join(OUT,"figures_pub","Fig_niche_model_summary.pdf"),bbox_inches="tight")
print("saved")
