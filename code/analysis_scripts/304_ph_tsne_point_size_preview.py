# -*- coding: utf-8 -*-
"""Preview-only PH t-SNE point-size comparison; no final export."""

from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.lines import Line2D

ROOT = Path(r"D:\桌面\single cell")
IN = ROOT / "analysis_output" / "figures_pub" / "2026-10-23_Fig2D_PH_tSNE_params" / "PH_tsne_v3_p50_pc20.csv"
OUT = ROOT / "analysis_output" / "figures_pub" / "2026-10-23_Fig2D_PH_tSNE_v3" / "preview_only"
OUT.mkdir(parents=True, exist_ok=True)

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
    "savefig.facecolor": "white",
})
colors = {0: "#4C78A8", 1: "#F58518", 2: "#E45756", 3: "#54A24B"}
names = {0: "s0 ECM/ribosome-high", 1: "s1 progenitor-like", 2: "s2 immune/ECM-primed", 3: "s3 TGM1-low"}
d = pd.read_csv(IN)
xmin, xmax = d["tsne1"].quantile([0.005, 0.995])
ymin, ymax = d["tsne2"].quantile([0.005, 0.995])
padx = (xmax - xmin) * 0.10
pady = (ymax - ymin) * 0.10
xlim = (xmin - padx, xmax + padx)
ylim = (ymin - pady, ymax + pady)

fig, axes = plt.subplots(1, 3, figsize=(8.2, 3.2))
for ax, size in zip(axes, [6, 9, 12]):
    for s in [0, 1, 2, 3]:
        z = d[d["subcluster"] == s]
        ax.scatter(z["tsne1"], z["tsne2"], s=size, alpha=0.72, color=colors[s], edgecolors="none", rasterized=True)
    for s in [0, 1, 2, 3]:
        z = d[d["subcluster"] == s]
        ax.text(z["tsne1"].median(), z["tsne2"].median(), f"s{s}", ha="center", va="center",
                fontsize=7.5, fontweight="bold", color="#1e293b",
                bbox=dict(boxstyle="round,pad=0.16", fc="white", ec=colors[s], lw=0.7, alpha=0.9))
    ax.set_xlim(xlim); ax.set_ylim(ylim)
    ax.set_box_aspect(1)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"point size = {size}", pad=4)
    for spine in ("top", "right"): ax.spines[spine].set_visible(False)

handles = [Line2D([0], [0], marker="o", linestyle="", markersize=5, markerfacecolor=colors[s], markeredgewidth=0, label=names[s]) for s in [0,1,2,3]]
fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.01), ncol=2, frameon=False, fontsize=6.4, handletextpad=0.4, columnspacing=1.2)
fig.tight_layout(rect=[0, 0.09, 1, 1])
fig.savefig(OUT / "PH_tsne_point_size_comparison_preview.png", dpi=300, bbox_inches="tight")
plt.close(fig)
print(OUT / "PH_tsne_point_size_comparison_preview.png")
