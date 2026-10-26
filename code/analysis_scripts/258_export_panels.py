# -*- coding: utf-8 -*-
"""258_export_panels.py
Export the six panels of the niche-model composite as standalone SVG and
1200-dpi TIFF files, using Times New Roman throughout.
"""
import os, sys
import numpy as np, pandas as pd, anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

OUT = r"D:\桌面\single cell\analysis_output"
S1 = os.path.join(OUT,"step1"); S3 = os.path.join(OUT,"step3")
DEST = os.path.join(OUT,"figures_pub","Fig6_niche_panels")
for sub in ["svg","tiff","png"]: os.makedirs(os.path.join(DEST,sub), exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", os.path.join(OUT,".mplconfig"))

FONT = "Times New Roman"
assert fm.findfont(FONT, fallback_to_default=False)
plt.rcParams.update({
    "font.family": "serif", "font.serif": [FONT], "font.size": 8,
    "axes.titlesize": 9, "axes.labelsize": 8, "xtick.labelsize": 7,
    "ytick.labelsize": 7, "legend.fontsize": 6.5,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
    "savefig.facecolor": "white",
})

a = ad.read_h5ad(os.path.join(S3,"step3_lyo_scores.h5ad"))
T = a.obsm["X_tsne"]
ov = pd.read_csv(os.path.join(S3,"step3_hypergeom_overlaps.csv"))
lr = pd.read_csv(os.path.join(S1,"step1_lr_communication_sonicparanoid.csv"))
order = ["PH","GH","TS","MH2","MyoEpi","M","F","R","B","PM-like"]
colors = {"PH":"#E45756","GH":"#72B7B2","TS":"#4C78A8","MH2":"#F58518","MyoEpi":"#54A24B",
          "M":"#B279A2","F":"#FF9DA6","R":"#9D755D","B":"#BAB0AC","PM-like":"#79706E"}

def save(fig, name):
    fig.savefig(os.path.join(DEST,"svg",f"{name}.svg"), format="svg", bbox_inches="tight")
    fig.savefig(os.path.join(DEST,"tiff",f"{name}.tiff"), format="tiff", dpi=1200,
                bbox_inches="tight", pil_kwargs={"compression":"tiff_lzw"})
    fig.savefig(os.path.join(DEST,"png",f"{name}.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("saved", name, flush=True)

def letter(ax, s):
    ax.text(-0.08, 1.04, s, transform=ax.transAxes, fontsize=10, fontweight="bold",
            va="top", ha="left")

# a cell types
fig, ax = plt.subplots(figsize=(3.6,3.3))
for ct in order:
    m = (a.obs.cell_type==ct).values
    ax.scatter(T[m,0],T[m,1],s=1.0,c=colors[ct],linewidths=0,rasterized=True)
    if m.sum()>20:
        ax.text(np.median(T[m,0]),np.median(T[m,1]),ct,fontsize=6,ha="center",
                bbox=dict(boxstyle="round,pad=0.10",fc="white",ec="none",alpha=0.65))
ax.set_xticks([]); ax.set_yticks([])
for s in ax.spines.values(): s.set_visible(False)
ax.set_xlabel("t-SNE 1", labelpad=2); ax.set_ylabel("t-SNE 2", labelpad=2)
letter(ax,"a"); fig.tight_layout(); save(fig,"panel_a_cell_types")

# b progenitor UCell t-SNE
fig, ax = plt.subplots(figsize=(3.9,3.3))
v = a.obs["ucell_LyO_Progenitor"].values
sc = ax.scatter(T[:,0],T[:,1],c=v,cmap="viridis",s=1.0,linewidths=0,rasterized=True,
                vmin=np.quantile(v,0.02),vmax=np.quantile(v,0.98))
ax.set_xticks([]); ax.set_yticks([])
for s in ax.spines.values(): s.set_visible(False)
ax.set_xlabel("t-SNE 1", labelpad=2); ax.set_ylabel("t-SNE 2", labelpad=2)
cb = plt.colorbar(sc, ax=ax, fraction=0.045, pad=0.01); cb.ax.tick_params(labelsize=6)
cb.set_label("UCell score", fontsize=7)
letter(ax,"b"); fig.tight_layout(); save(fig,"panel_b_progenitor_tsne")

# c progenitor boxplot
fig, ax = plt.subplots(figsize=(3.8,3.0))
data = [a.obs.loc[a.obs.cell_type==c,"ucell_LyO_Progenitor"].values for c in order]
bp = ax.boxplot(data, tick_labels=order, showfliers=False, patch_artist=True, widths=0.62,
                medianprops=dict(color="#222222", lw=0.7), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
for patch,c in zip(bp["boxes"],order): patch.set_facecolor(colors[c]); patch.set_alpha(0.75); patch.set_edgecolor("#333333"); patch.set_linewidth(0.5)
ax.set_ylabel("UCell score"); ax.tick_params(axis="x", rotation=60)
for s in ["top","right"]: ax.spines[s].set_visible(False)
letter(ax,"c"); fig.tight_layout(); save(fig,"panel_c_progenitor_boxplot")

# d overlap bar
fig, ax = plt.subplots(figsize=(4.0,2.8))
sub = ov[ov.lyo_set=="LyO_Progenitor"].set_index("compared_to").reindex(
    ["PH_markers","niche_markers","TS_markers","MH2_markers","MyoEpi_markers","M_markers"])
y = np.arange(len(sub)); vals = sub["fold"].fillna(0).values
ax.barh(y, vals, color=["#C44E52" if x>1 else "#BBBBBB" for x in vals], height=0.65)
ax.set_yticks(y); ax.set_yticklabels(["PH","niche","TS","MH2","MyoEpi","M"]); ax.invert_yaxis()
for i,(val,p,k) in enumerate(zip(vals,sub["p_hypergeom"].values,sub["overlap"].values)):
    ax.text(max(val,0.02)+2, i, f"{int(k)} genes, {val:.1f}×\nP = {p:.1e}", va="center", fontsize=6)
ax.set_xlim(0, max(12, vals.max()*1.7)); ax.set_xlabel("Fold enrichment")
for s in ["top","right"]: ax.spines[s].set_visible(False)
letter(ax,"d"); fig.tight_layout(); save(fig,"panel_d_progenitor_overlap")

# e LR heatmap
fig, ax = plt.subplots(figsize=(3.1,2.8))
piv = lr.groupby(["source","target"])["prob"].sum().unstack(fill_value=0).reindex(
    index=["TS","MH2","MyoEpi","M"], columns=["PH","GH"], fill_value=0)
im = ax.imshow(piv.values, cmap="Reds", aspect="auto")
ax.set_xticks([0,1]); ax.set_xticklabels(["PH","GH"]); ax.set_yticks(range(4)); ax.set_yticklabels(["TS","MH2","MyoEpi","M"])
for i in range(4):
    for j in range(2):
        ax.text(j,i,f"{piv.values[i,j]:.2f}",ha="center",va="center",fontsize=7,
                color="white" if piv.values[i,j]>piv.values.max()*0.6 else "#333333")
cb = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.02); cb.ax.tick_params(labelsize=6)
cb.set_label("Summed probability", fontsize=7)
letter(ax,"e"); fig.tight_layout(); save(fig,"panel_e_lr_heatmap")

# f projected pseudotime
fig, ax = plt.subplots(figsize=(4.0,3.0))
data = [a.obs.loc[a.obs.cell_type==c,"projected_pseudotime"].values for c in order]
bp = ax.boxplot(data, tick_labels=order, showfliers=False, patch_artist=True, widths=0.62,
                medianprops=dict(color="#222222", lw=0.7), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
for patch,c in zip(bp["boxes"],order): patch.set_facecolor(colors[c]); patch.set_alpha(0.75); patch.set_edgecolor("#333333"); patch.set_linewidth(0.5)
ax.set_ylabel("Projected pseudotime"); ax.tick_params(axis="x", rotation=60)
for s in ["top","right"]: ax.spines[s].set_visible(False)
letter(ax,"f"); fig.tight_layout(); save(fig,"panel_f_projected_pseudotime")

print("all panels exported to", DEST)
