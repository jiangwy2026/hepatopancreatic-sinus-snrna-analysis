# -*- coding: utf-8 -*-
"""Preview a PH-only t-SNE panel with substate composition shift."""

from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.lines import Line2D

ROOT = Path(r"D:\桌面\single cell")
OUT = ROOT / "analysis_output" / "figures_pub" / "2026-10-23_Fig2D_PH_tSNE"
OUT.mkdir(parents=True, exist_ok=True)
H5 = ROOT / "analysis_output" / "step3" / "step3_lyo_scores.h5ad"
SUB = ROOT / "analysis_output" / "dynamic_subpop_subclusters" / "PH_subcluster_umap.csv"
COMP = ROOT / "analysis_output" / "Fig5_dynamic_subpop_subcluster_composition.csv"

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
names = {
    0: "s0 ECM/ribosome-high",
    1: "s1 progenitor-like",
    2: "s2 immune/ECM-primed",
    3: "s3 TGM1-low",
}
order = [0, 1, 2, 3]

a = ad.read_h5ad(H5, backed="r")
meta = a.obs.loc[a.obs["cell_type"].astype(str) == "PH", ["group"]].copy()
tsne = np.asarray(a.obsm["X_tsne"])[a.obs["cell_type"].astype(str).eq("PH").to_numpy(), :]
a.file.close()
meta["cell"] = meta.index.astype(str)
meta["group"] = meta["group"].astype(str).replace({"blank": "Control", "treated": "Hemospasia"})

sub = pd.read_csv(SUB, index_col=0)
sub.index = sub.index.astype(str)
sub = sub[["subcluster"]]
d = meta.merge(sub, left_on="cell", right_index=True, how="inner")
coords = pd.DataFrame(tsne, index=meta.index.astype(str), columns=["tsne1", "tsne2"])
d = d.merge(coords, left_on="cell", right_index=True, how="inner")
d["subcluster"] = d["subcluster"].astype(int)
print("PH cells mapped", len(d), "groups", d["group"].value_counts().to_dict(), "subclusters", d["subcluster"].value_counts().to_dict())

comp = pd.read_csv(COMP)
comp = comp[comp.population == "PH"].sort_values("subcluster").set_index("subcluster")

fig = plt.figure(figsize=(7.7, 3.25))
gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 0.95], wspace=0.32, left=0.06, right=0.97, bottom=0.18, top=0.86)

# Shared t-SNE limits
xmin, xmax = d["tsne1"].min() - 0.5, d["tsne1"].max() + 0.5
ymin, ymax = d["tsne2"].min() - 0.5, d["tsne2"].max() + 0.5

for j, grp in enumerate(["Control", "Hemospasia"]):
    ax = fig.add_subplot(gs[0, j])
    subd = d[d["group"] == grp]
    for s in order:
        ss = subd[subd["subcluster"] == s]
        ax.scatter(ss["tsne1"], ss["tsne2"], s=5.0, color=colors[s], alpha=0.72,
                   edgecolors="none", rasterized=True)
    ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_box_aspect(1)
    ax.set_xlabel("t-SNE 1")
    if j == 0: ax.set_ylabel("t-SNE 2")
    ax.set_title(f"{grp} (n={len(subd)})", pad=5)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

# Legend under t-SNE panels
handles = [Line2D([0], [0], marker="o", linestyle="", markersize=5,
                  markerfacecolor=colors[s], markeredgewidth=0, label=names[s]) for s in order]
fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.07, 0.015),
           ncol=2, frameon=False, fontsize=6.4, handletextpad=0.4, columnspacing=1.2)

ax = fig.add_subplot(gs[0, 2])
for s in order:
    c = float(comp.loc[s, "pct_control"])
    h = float(comp.loc[s, "pct_hemospasia"])
    fold = float(comp.loc[s, "fold_hemospasia_vs_control"])
    ax.plot([c, h], [s, s], color="#BDBDBD", lw=1.6, zorder=1)
    ax.scatter([c], [s], s=34, color="#8A8A8A", zorder=3)
    ax.scatter([h], [s], s=34, color=colors[s], zorder=3)
    ax.text(105, s, f"{fold:.2f}×", ha="left", va="center", fontsize=6.6, color=colors[s], fontweight="bold")
ax.set_yticks(order)
ax.set_yticklabels([f"s{s}" for s in order])
ax.set_xlim(0, 125)
ax.set_xticks([0, 20, 40, 60, 80, 100])
ax.set_xlabel("PH cells (%)")
ax.set_title("Substate shift", pad=22)
ax.legend(handles=[Line2D([0],[0], marker="o", linestyle="", markersize=5, markerfacecolor="#8A8A8A", markeredgewidth=0, label="Control"),
                   Line2D([0],[0], marker="o", linestyle="", markersize=5, markerfacecolor="#E45756", markeredgewidth=0, label="Hemospasia")],
          frameon=False, fontsize=6.1, loc="lower center", bbox_to_anchor=(0.45, 1.01), ncol=2)
for spine in ("top", "right"):
    ax.spines[spine].set_visible(False)

fig.savefig(OUT / "Fig2D_PH_tSNE_preview.png", dpi=300, bbox_inches="tight")
fig.savefig(OUT / "Fig2D_PH_tSNE_preview.svg", format="svg", bbox_inches="tight")
fig.savefig(OUT / "Fig2D_PH_tSNE_preview.tiff", format="tiff", dpi=1200, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
plt.close(fig)
print(OUT / "Fig2D_PH_tSNE_preview.png")
