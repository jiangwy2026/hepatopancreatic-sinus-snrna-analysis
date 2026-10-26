# -*- coding: utf-8 -*-
"""Compute and preview a PH-only t-SNE using the Harmony PCA."""

from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.lines import Line2D

ROOT = Path(r"D:\桌面\single cell")
OUT = ROOT / "analysis_output" / "figures_pub" / "2026-10-23_Fig2D_PH_tSNE_v2"
OUT.mkdir(parents=True, exist_ok=True)
H5 = ROOT / "analysis_output" / "step3" / "step3_lyo_scores.h5ad"
SUB = ROOT / "analysis_output" / "dynamic_subpop_subclusters" / "PH_subcluster_umap.csv"

assert fm.findfont("Times New Roman", fallback_to_default=False)
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman"],
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "axes.linewidth": 0.65,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "savefig.facecolor": "white",
})

colors = {0: "#4C78A8", 1: "#F58518", 2: "#E45756", 3: "#54A24B"}
names = {0: "s0 ECM/ribosome-high", 1: "s1 progenitor-like", 2: "s2 immune/ECM-primed", 3: "s3 TGM1-low"}

a = ad.read_h5ad(H5)
ph = a.obs["cell_type"].astype(str).eq("PH").to_numpy()
a = a[ph].copy()
sub = pd.read_csv(SUB, index_col=0)
sub.index = sub.index.astype(str)
a = a[a.obs_names.astype(str)]
a.obs["subcluster"] = sub.reindex(a.obs_names.astype(str))["subcluster"].astype(int).values
a.obs["group_label"] = a.obs["group"].astype(str).replace({"blank": "Control", "treated": "Hemospasia"})
print("PH shape", a.shape, "subclusters", a.obs["subcluster"].value_counts().sort_index().to_dict())

# The Harmony PCA was computed on all Hs cells; recompute the neighborhood and t-SNE only within PH.
sc.pp.neighbors(a, use_rep="X_pca_harmony", n_neighbors=20, metric="cosine", random_state=0)
sc.tl.tsne(
    a, use_rep="X_pca_harmony", n_pcs=None, perplexity=30,
    learning_rate="auto", early_exaggeration=12, random_state=42, n_jobs=4,
)
a.write_h5ad(OUT / "PH_only_tsne.h5ad", compression="gzip")

xy = pd.DataFrame(a.obsm["X_tsne"], columns=["tsne1", "tsne2"], index=a.obs_names)
xy["group"] = a.obs["group_label"].astype(str).values
xy["subcluster"] = a.obs["subcluster"].astype(int).values
xy.to_csv(OUT / "PH_only_tsne_coordinates.csv", encoding="utf-8-sig")

xmin, xmax = xy["tsne1"].quantile([0.005, 0.995])
ymin, ymax = xy["tsne2"].quantile([0.005, 0.995])
padx = (xmax - xmin) * 0.10
pady = (ymax - ymin) * 0.10
xlim = (xmin - padx, xmax + padx)
ylim = (ymin - pady, ymax + pady)

fig = plt.figure(figsize=(8.2, 4.3))
gs = fig.add_gridspec(1, 2, width_ratios=[1.35, 1.0], wspace=0.28, left=0.06, right=0.98, bottom=0.18, top=0.88)

# Combined PH t-SNE
ax = fig.add_subplot(gs[0, 0])
for s in [0, 1, 2, 3]:
    z = xy[xy["subcluster"] == s]
    ax.scatter(z["tsne1"], z["tsne2"], s=5.0, alpha=0.72, color=colors[s], edgecolors="none", rasterized=True)
for s in [0, 1, 2, 3]:
    z = xy[xy["subcluster"] == s]
    cx, cy = z["tsne1"].median(), z["tsne2"].median()
    ax.text(cx, cy, f"s{s}", ha="center", va="center", fontsize=8, fontweight="bold",
            color="#1e293b", bbox=dict(boxstyle="round,pad=0.18", fc="white", ec=colors[s], lw=0.8, alpha=0.9))
ax.set_xlim(xlim); ax.set_ylim(ylim)
ax.set_xlabel("PH-only t-SNE 1"); ax.set_ylabel("PH-only t-SNE 2")
ax.set_title("PH substate structure", pad=5)
ax.set_box_aspect(1)
for spine in ("top", "right"): ax.spines[spine].set_visible(False)

# Control/Hemospasia split
right = fig.add_gridspec(2, 1, left=0.55, right=0.98, bottom=0.18, top=0.88, hspace=0.32)
for j, grp in enumerate(["Control", "Hemospasia"]):
    ax = fig.add_subplot(right[j, 0])
    z = xy[xy["group"] == grp]
    for s in [0, 1, 2, 3]:
        zz = z[z["subcluster"] == s]
        ax.scatter(zz["tsne1"], zz["tsne2"], s=4.0, alpha=0.68, color=colors[s], edgecolors="none", rasterized=True)
    ax.set_xlim(xlim); ax.set_ylim(ylim)
    ax.set_box_aspect(1)
    ax.set_title(f"{grp} (n={len(z)})", pad=3)
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ("top", "right"): ax.spines[spine].set_visible(False)

handles = [Line2D([0], [0], marker="o", linestyle="", markersize=5, markerfacecolor=colors[s], markeredgewidth=0, label=names[s]) for s in [0,1,2,3]]
fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.06, 0.02), ncol=2, frameon=False, fontsize=6.2, handletextpad=0.4, columnspacing=1.1)
fig.text(0.55, 0.91, "Group overlay on the same PH-only embedding", ha="center", va="bottom", fontsize=9, fontweight="bold")

fig.savefig(OUT / "Fig2D_PH_only_tSNE_preview.png", dpi=300, bbox_inches="tight")
fig.savefig(OUT / "Fig2D_PH_only_tSNE_preview.svg", format="svg", bbox_inches="tight")
fig.savefig(OUT / "Fig2D_PH_only_tSNE_preview.tiff", format="tiff", dpi=1200, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
plt.close(fig)
print(OUT / "Fig2D_PH_only_tSNE_preview.png")
