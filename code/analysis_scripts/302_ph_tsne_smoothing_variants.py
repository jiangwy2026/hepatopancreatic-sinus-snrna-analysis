# -*- coding: utf-8 -*-
"""Compare PH-only t-SNE smoothing parameter sets."""

from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

ROOT = Path(r"D:\桌面\single cell")
OUT = ROOT / "analysis_output" / "figures_pub" / "2026-10-23_Fig2D_PH_tSNE_params"
OUT.mkdir(parents=True, exist_ok=True)
BASE = ROOT / "analysis_output" / "figures_pub" / "2026-10-23_Fig2D_PH_tSNE_v2" / "PH_only_tsne.h5ad"

assert fm.findfont("Times New Roman", fallback_to_default=False)
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman"],
    "font.size": 7.5,
    "axes.titlesize": 8.5,
    "axes.labelsize": 7.5,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "axes.linewidth": 0.6,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "savefig.facecolor": "white",
})
colors = {0: "#4C78A8", 1: "#F58518", 2: "#E45756", 3: "#54A24B"}
names = {0: "s0 ECM/ribosome-high", 1: "s1 progenitor-like", 2: "s2 immune/ECM-primed", 3: "s3 TGM1-low"}

variants = [
    ("v1_p40_pc30", 40, 30, 12),
    ("v2_p60_pc30", 60, 30, 12),
    ("v3_p50_pc20", 50, 20, 8),
    ("v4_p80_pc20", 80, 20, 8),
]

a0 = ad.read_h5ad(BASE)
all_xy = {}
for name, perp, npc, ee in variants:
    a = a0.copy()
    print("running", name, "perplexity", perp, "n_pcs", npc, "ee", ee, flush=True)
    sc.tl.tsne(
        a, use_rep="X_pca_harmony", n_pcs=npc, perplexity=perp,
        learning_rate="auto", early_exaggeration=ee, random_state=42, n_jobs=4,
    )
    xy = pd.DataFrame(a.obsm["X_tsne"], columns=["tsne1", "tsne2"], index=a.obs_names)
    xy["subcluster"] = a.obs["subcluster"].astype(int).values
    xy["group"] = a.obs["group_label"].astype(str).values
    all_xy[name] = xy
    xy.to_csv(OUT / f"PH_tsne_{name}.csv", encoding="utf-8-sig")

# Plot comparison grid
fig, axes = plt.subplots(2, 2, figsize=(7.1, 6.2))
for ax, (name, perp, npc, ee) in zip(axes.ravel(), variants):
    xy = all_xy[name]
    for s in [0, 1, 2, 3]:
        z = xy[xy["subcluster"] == s]
        ax.scatter(z["tsne1"], z["tsne2"], s=4.0, alpha=0.72, color=colors[s], edgecolors="none", rasterized=True)
    for s in [0, 1, 2, 3]:
        z = xy[xy["subcluster"] == s]
        ax.text(z["tsne1"].median(), z["tsne2"].median(), f"s{s}", ha="center", va="center",
                fontsize=6.5, fontweight="bold", color="#1e293b",
                bbox=dict(boxstyle="round,pad=0.16", fc="white", ec=colors[s], lw=0.6, alpha=0.9))
    ax.set_title(f"{name}: perplexity {perp}, PCs {npc}, EE {ee}", pad=4)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_box_aspect(1)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)

handles = [plt.Line2D([0], [0], marker="o", linestyle="", markersize=5, markerfacecolor=colors[s], markeredgewidth=0, label=names[s]) for s in [0,1,2,3]]
fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.01), ncol=2, frameon=False, fontsize=6.5, handletextpad=0.4, columnspacing=1.2)
fig.tight_layout(rect=[0, 0.06, 1, 1])
fig.savefig(OUT / "PH_tsne_parameter_comparison.png", dpi=300, bbox_inches="tight")
fig.savefig(OUT / "PH_tsne_parameter_comparison.svg", format="svg", bbox_inches="tight")
plt.close(fig)
print(OUT / "PH_tsne_parameter_comparison.png")
