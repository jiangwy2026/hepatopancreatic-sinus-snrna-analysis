# -*- coding: utf-8 -*-
"""Supplementary Fig. 14: marker-level LyO progenitor program mirror.

This analysis uses the available LyO top50 markers only. It does not claim
cell-level UCell scores or cell-cycle adjustment in the LyO compartment.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import Rectangle

ROOT = Path(r"D:\桌面\single cell")
OUT = ROOT / "analysis_output"
DEST = OUT / "figures_pub" / "Supplementary_niche"
for sub in ("png", "svg", "tiff", "pdf"):
    (DEST / sub).mkdir(parents=True, exist_ok=True)

assert fm.findfont("Times New Roman", fallback_to_default=False)
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman"],
    "font.size": 7,
    "axes.titlesize": 8,
    "axes.labelsize": 7,
    "xtick.labelsize": 6.2,
    "ytick.labelsize": 6.2,
    "axes.linewidth": 0.6,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "savefig.facecolor": "white",
})

# Available data: LyO top50 marker table, cross-organ marker mapping and Hs UCell scores.
heat = pd.read_csv(OUT / "lyo_cross_organ_mapping_scores.csv")
heat = heat[heat["organ"] == "Hepatic sinus"].copy()
best = pd.read_csv(OUT / "Fig_lyo_cluster_best_matches.csv")
hyper = pd.read_csv(OUT / "step3" / "step3_hypergeom_overlaps.csv")

cell_order = ["PH", "GH", "TS", "MH2", "F", "R", "B", "M", "MyoEpi", "PM-like"]
name_map = dict(zip(best.lyo_cluster, best.lyo_name))
mat = heat.pivot(index="lyo_cluster", columns="target_group", values="z_score").reindex(columns=cell_order)
mat = mat.sort_index()
cluster_labels = [f"{int(c)} {name_map.get(int(c), '')}" for c in mat.index]

# Hs per-cell LyO progenitor UCell score
a = ad.read_h5ad(OUT / "step3" / "step3_lyo_scores.h5ad", backed="r")
obs = a.obs.copy()
a.file.close()
ucell = obs.groupby("cell_type")["ucell_LyO_Progenitor"].agg(["count", "mean", "median"])
ucell["frac_above_global"] = obs.groupby("cell_type")["ucell_LyO_Progenitor"].apply(
    lambda x: float((x > obs["ucell_LyO_Progenitor"].median()).mean())
)
ucell = ucell.reindex(cell_order).dropna(subset=["mean"])

# Panel C: progenitor program overlaps with Hs marker programs
sub = hyper[hyper["lyo_set"] == "LyO_Progenitor"].copy()
sub["label"] = sub["compared_to"].str.replace("_markers", "", regex=False).replace({"niche": "Niche"})
sub = sub.set_index("label").reindex(["PH", "Niche", "TS", "MH2", "MyoEpi", "M"]).reset_index()

fig = plt.figure(figsize=(7.4, 5.0))
gs = fig.add_gridspec(2, 2, width_ratios=[1.65, 1.0], height_ratios=[1.0, 0.82], wspace=0.38, hspace=0.48,
                      left=0.08, right=0.97, top=0.92, bottom=0.10)

# A: marker-level mapping heatmap
ax = fig.add_subplot(gs[0, 0])
vmax = np.nanmax(np.abs(mat.values))
im = ax.imshow(mat.values, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
ax.set_xticks(range(len(cell_order))); ax.set_xticklabels(cell_order, rotation=55, ha="right")
ax.set_yticks(range(len(cluster_labels))); ax.set_yticklabels(cluster_labels, fontsize=5.4)
ax.set_xlabel("Hs cell state")
ax.set_ylabel("Lymphoid-organ cluster")
ax.set_title("LyO cluster-to-Hs marker mapping", pad=4)
# highlight progenitor cluster 13
if 13 in mat.index:
    row = list(mat.index).index(13)
    ax.add_patch(Rectangle((-0.5, row - 0.5), len(cell_order), 1, fill=False, ec="#111111", lw=1.05))
for s in ax.spines.values(): s.set_visible(False)
ax.tick_params(length=0, pad=2)
ax.text(-0.12, 1.03, "A", transform=ax.transAxes, fontsize=10, fontweight="bold")
cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
cb.ax.tick_params(labelsize=5.4, length=1.5)
cb.set_label("marker-overlap z", fontsize=5.8)

# B: LyO progenitor UCell score across Hs states
ax = fig.add_subplot(gs[0, 1])
y = np.arange(len(ucell))
ax.hlines(y, 0.46, ucell["mean"].values, color="#D0D0D0", lw=1.0)
ax.scatter(ucell["mean"].values, y, s=20 + 40 * ucell["frac_above_global"].values,
           color=["#B2182B" if x == "PH" else "#6B6B6B" for x in ucell.index], zorder=3)
ax.set_yticks(y); ax.set_yticklabels(ucell.index)
ax.invert_yaxis()
ax.set_xlim(0.46, 0.555)
ax.set_xlabel("Mean LyO progenitor UCell")
ax.set_title("LyO progenitor program in Hs", pad=4)
for i, (name, r) in enumerate(ucell.iterrows()):
    ax.text(0.531, i, f"{r['mean']:.3f} / {100*r['frac_above_global']:.0f}%",
            va="center", ha="left", fontsize=4.9, color="#B2182B" if name == "PH" else "#444444")
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.text(-0.18, 1.03, "B", transform=ax.transAxes, fontsize=10, fontweight="bold")

# C: hypergeometric overlap of LyO progenitor top50 with Hs programs
ax = fig.add_subplot(gs[1, :])
y = np.arange(len(sub))[::-1]
cols = ["#B2182B" if x == "PH" else "#8A8A8A" for x in sub["label"]]
ax.hlines(y, 0, sub["fold"].values, color="#D0D0D0", lw=1.0)
ax.scatter(sub["fold"].values, y, s=35, color=cols, zorder=3)
for yi, (_, r) in zip(y, sub.iterrows()):
    p = r["p_hypergeom"]
    pstr = "P=1" if p >= 0.999 else f"P={p:.1e}"
    ax.text(r["fold"] + 1.8, yi, f"{int(r['overlap'])} genes; {r['fold']:.1f}×; {pstr}",
            va="center", fontsize=6.0, color="#B2182B" if r["label"] == "PH" else "#444444")
ax.set_yticks(y); ax.set_yticklabels(sub["label"])
ax.set_xlim(0, 80)
ax.set_xlabel("Fold enrichment of LyO progenitor top50 in Hs marker programs")
ax.set_title("Marker-level program overlap", pad=4)
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.text(-0.06, 1.03, "C", transform=ax.transAxes, fontsize=10, fontweight="bold")
ax.text(0.99, -0.32, "Marker-level only; full LyO matrix unavailable; no cell-cycle adjustment",
        transform=ax.transAxes, ha="right", va="top", fontsize=5.6, color="#666666")

# Source data
mat.to_csv(DEST / "Supp14_LyO_Hs_mapping_matrix.csv", encoding="utf-8-sig")
ucell.to_csv(DEST / "Supp14_LyO_progenitor_Hs_UCell.csv", encoding="utf-8-sig")
sub.to_csv(DEST / "Supp14_LyO_progenitor_overlap.csv", index=False, encoding="utf-8-sig")

stem = DEST / "Supp14_LyO_marker_mirror"
fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
fig.savefig(stem.with_suffix(".svg"), format="svg", bbox_inches="tight")
fig.savefig(stem.with_suffix(".pdf"), format="pdf", bbox_inches="tight")
fig.savefig(stem.with_suffix(".tiff"), format="tiff", dpi=1200, bbox_inches="tight",
            pil_kwargs={"compression": "tiff_lzw"})
plt.close(fig)
print(stem)
