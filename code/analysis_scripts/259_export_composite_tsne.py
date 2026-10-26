# -*- coding: utf-8 -*-
"""259_export_composite_tsne.py — full 6-panel composite, Times New Roman, SVG + 1200 dpi TIFF."""
import os
import numpy as np, pandas as pd, anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

OUT=r"D:\桌面\single cell\analysis_output"; S1=os.path.join(OUT,"step1"); S3=os.path.join(OUT,"step3")
DEST=os.path.join(OUT,"figures_pub","Fig6_niche_panels")
os.environ.setdefault("MPLCONFIGDIR",os.path.join(OUT,".mplconfig"))
assert fm.findfont("Times New Roman", fallback_to_default=False)
plt.rcParams.update({"font.family":"serif","font.serif":["Times New Roman"],"font.size":7,
 "axes.titlesize":8,"axes.labelsize":7,"xtick.labelsize":6.5,"ytick.labelsize":6.5,
 "axes.linewidth":0.6,"svg.fonttype":"none","savefig.facecolor":"white"})

a=ad.read_h5ad(os.path.join(S3,"step3_lyo_scores.h5ad")); T=a.obsm["X_tsne"]
ov=pd.read_csv(os.path.join(S3,"step3_hypergeom_overlaps.csv"))
lr=pd.read_csv(os.path.join(S1,"step1_lr_communication_sonicparanoid.csv"))
order=["PH","GH","TS","MH2","MyoEpi","M","F","R","B","PM-like"]
colors={"PH":"#E45756","GH":"#72B7B2","TS":"#4C78A8","MH2":"#F58518","MyoEpi":"#54A24B",
 "M":"#B279A2","F":"#FF9DA6","R":"#9D755D","B":"#BAB0AC","PM-like":"#79706E"}
fig=plt.figure(figsize=(7.2,7.0)); gs=fig.add_gridspec(3,2,hspace=0.42,wspace=0.28,left=0.08,right=0.97,top=0.95,bottom=0.07)
ax=fig.add_subplot(gs[0,0])
for ct in order:
    m=(a.obs.cell_type==ct).values
    ax.scatter(T[m,0],T[m,1],s=1.2,c=colors[ct],linewidths=0,rasterized=True)
    if m.sum()>20: ax.text(np.median(T[m,0]),np.median(T[m,1]),ct,fontsize=6,ha="center",bbox=dict(boxstyle="round,pad=0.12",fc="white",ec="none",alpha=0.65))
ax.set_title("a  Cell types (t-SNE)"); ax.set_xticks([]); ax.set_yticks([])
for s in ax.spines.values(): s.set_visible(False)
ax=fig.add_subplot(gs[0,1]); v=a.obs["ucell_LyO_Progenitor"].values
sc=ax.scatter(T[:,0],T[:,1],c=v,cmap="viridis",s=1.2,linewidths=0,rasterized=True,vmin=np.quantile(v,0.02),vmax=np.quantile(v,0.98))
ax.set_title("b  LyO progenitor program (UCell, t-SNE)"); ax.set_xticks([]); ax.set_yticks([])
for s in ax.spines.values(): s.set_visible(False)
plt.colorbar(sc,ax=ax,fraction=0.042,pad=0.01)
ax=fig.add_subplot(gs[1,0]); data=[a.obs.loc[a.obs.cell_type==c,"ucell_LyO_Progenitor"].values for c in order]
bp=ax.boxplot(data,tick_labels=order,showfliers=False,patch_artist=True,widths=0.62,medianprops=dict(color="#222",lw=0.7))
for p,c in zip(bp["boxes"],order): p.set_facecolor(colors[c]); p.set_alpha(0.75); p.set_edgecolor("#333"); p.set_linewidth(0.5)
ax.set_title("c  Progenitor program by cell type"); ax.tick_params(axis="x",rotation=60); ax.set_ylabel("UCell score")
ax=fig.add_subplot(gs[1,1])
sub=ov[ov.lyo_set=="LyO_Progenitor"].set_index("compared_to").reindex(["PH_markers","niche_markers","TS_markers","MH2_markers","MyoEpi_markers","M_markers"])
y=np.arange(len(sub)); vals=sub["fold"].fillna(0).values
ax.barh(y,vals,color=["#C44E52" if x>1 else "#BBB" for x in vals])
ax.set_yticks(y); ax.set_yticklabels(["PH","niche","TS","MH2","MyoEpi","M"]); ax.invert_yaxis()
for i,(val,p,k) in enumerate(zip(vals,sub["p_hypergeom"].values,sub["overlap"].values)):
    ax.text(max(val,0.02)+2,i,f"{int(k)} genes, {val:.1f}×\nP = {p:.1e}",va="center",fontsize=5.5)
ax.set_xlim(0,max(12,vals.max()*1.6)); ax.set_xlabel("Fold enrichment"); ax.set_title("d  Progenitor top50 overlap")
for s in ["top","right"]: ax.spines[s].set_visible(False)
ax=fig.add_subplot(gs[2,0])
piv=lr.groupby(["source","target"])["prob"].sum().unstack(fill_value=0).reindex(index=["TS","MH2","MyoEpi","M"],columns=["PH","GH"],fill_value=0)
im=ax.imshow(piv.values,cmap="Reds",aspect="auto"); ax.set_xticks([0,1]); ax.set_xticklabels(["PH","GH"])
ax.set_yticks(range(4)); ax.set_yticklabels(["TS","MH2","MyoEpi","M"])
for i in range(4):
    for j in range(2):
        ax.text(j,i,f"{piv.values[i,j]:.2f}",ha="center",va="center",fontsize=6,color="white" if piv.values[i,j]>piv.values.max()*0.6 else "#333")
ax.set_title("e  Niche -> hemocyte ECM-integrin"); plt.colorbar(im,ax=ax,fraction=0.042,pad=0.01)
ax=fig.add_subplot(gs[2,1]); data=[a.obs.loc[a.obs.cell_type==c,"projected_pseudotime"].values for c in order]
bp=ax.boxplot(data,tick_labels=order,showfliers=False,patch_artist=True,widths=0.62,medianprops=dict(color="#222",lw=0.7))
for p,c in zip(bp["boxes"],order): p.set_facecolor(colors[c]); p.set_alpha(0.75); p.set_edgecolor("#333"); p.set_linewidth(0.5)
ax.set_title("f  Projection onto PH->GH trajectory"); ax.tick_params(axis="x",rotation=60); ax.set_ylabel("Projected pseudotime")
fig.savefig(os.path.join(DEST,"svg","Fig6_niche_composite.svg"),format="svg",bbox_inches="tight")
fig.savefig(os.path.join(DEST,"tiff","Fig6_niche_composite.tiff"),format="tiff",dpi=1200,bbox_inches="tight",pil_kwargs={"compression":"tiff_lzw"})
fig.savefig(os.path.join(DEST,"png","Fig6_niche_composite.png"),dpi=200,bbox_inches="tight")
print("composite saved")
