# -*- coding: utf-8 -*-
"""Preview-only three equal PH-only t-SNE panels."""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch
from scipy.interpolate import splprep, splev

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
trajectory_order = [1, 0, 2, 3]
d = pd.read_csv(IN)
centroids = {s: (d.loc[d.subcluster == s, "tsne1"].median(), d.loc[d.subcluster == s, "tsne2"].median()) for s in [0,1,2,3]}
pts = np.array([centroids[s] for s in trajectory_order], dtype=float)
tck, _ = splprep([pts[:,0], pts[:,1]], s=0, k=2)
u = np.linspace(0, 1, 160)
sx, sy = splev(u, tck)
curve = np.column_stack([sx, sy])

xmin, xmax = d["tsne1"].quantile([0.005, 0.995])
ymin, ymax = d["tsne2"].quantile([0.005, 0.995])
padx = (xmax - xmin) * 0.10
pady = (ymax - ymin) * 0.10
xlim = (xmin - padx, xmax + padx)
ylim = (ymin - pady, ymax + pady)

fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.0))
groups = [("All PH", d), ("Control", d[d.group == "Control"]), ("Hemospasia", d[d.group == "Hemospasia"])]

for j, (title, data) in enumerate(groups):
    ax = axes[j]
    for s in [0, 1, 2, 3]:
        z = data[data.subcluster == s]
        ax.scatter(z["tsne1"], z["tsne2"], s=12, alpha=0.70, color=colors[s], edgecolors="none", rasterized=True)
    if j == 0:
        ax.plot(curve[:,0], curve[:,1], color="white", lw=4.0, alpha=0.95, zorder=5)
        ax.plot(curve[:,0], curve[:,1], color="#111111", lw=1.8, zorder=6)
        ax.add_patch(FancyArrowPatch(curve[-12], curve[-1], arrowstyle="-|>", mutation_scale=10,
                                     lw=1.2, color="#111111", zorder=7))
        for st in trajectory_order:
            cx, cy = centroids[st]
            ax.text(cx, cy, f"s{st}", ha="center", va="center", fontsize=7.2,
                    fontweight="bold", color="#1e293b", zorder=9,
                    bbox=dict(boxstyle="round,pad=0.16", fc="white", ec=colors[st],
                              lw=0.65, alpha=0.95))
    ax.set_xlim(xlim); ax.set_ylim(ylim)
    ax.set_box_aspect(1)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("tSNE")
    ax.set_ylabel("tSNE")
    ax.set_title(title, pad=4)
    for spine in ("top", "right"): ax.spines[spine].set_visible(False)

handles = [Line2D([0], [0], marker="o", linestyle="", markersize=5, markerfacecolor=colors[s], markeredgewidth=0, label=names[s]) for s in [0,1,2,3]]
fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.01), ncol=2, frameon=False, fontsize=6.3, handletextpad=0.4, columnspacing=1.1)
fig.tight_layout(rect=[0, 0.10, 1, 1])
fig.savefig(OUT / "PH_tsne_three_equal_panels_preview.png", dpi=300, bbox_inches="tight")
plt.close(fig)
print(OUT / "PH_tsne_three_equal_panels_preview.png")
