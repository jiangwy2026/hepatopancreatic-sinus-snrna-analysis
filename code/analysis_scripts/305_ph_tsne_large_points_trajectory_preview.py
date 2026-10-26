# -*- coding: utf-8 -*-
"""Preview-only PH t-SNE with larger points and s1->s0->s2->s3 centroid trajectory."""

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
# Smooth curve through the four centroids.
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

fig = plt.figure(figsize=(8.2, 4.3))
gs = fig.add_gridspec(1, 2, width_ratios=[1.35, 1.0], wspace=0.28, left=0.06, right=0.98, bottom=0.18, top=0.88)

def draw_trajectory(ax):
    ax.plot(curve[:,0], curve[:,1], color="white", lw=4.8, alpha=0.95, zorder=5)
    ax.plot(curve[:,0], curve[:,1], color="#111111", lw=2.1, zorder=6)
    # arrow at the end of the curve
    ax.add_patch(FancyArrowPatch(curve[-12], curve[-1], arrowstyle="-|>",
                                 mutation_scale=12, lw=1.5, color="#111111", zorder=7))

ax = fig.add_subplot(gs[0, 0])
for s in [0, 1, 2, 3]:
    z = d[d["subcluster"] == s]
    ax.scatter(z["tsne1"], z["tsne2"], s=14, alpha=0.68, color=colors[s], edgecolors="none", rasterized=True)
for s in [0, 1, 2, 3]:
    x, y = centroids[s]
    ax.text(x, y, f"s{s}", ha="center", va="center", fontsize=8, fontweight="bold",
            color="#1e293b", bbox=dict(boxstyle="round,pad=0.18", fc="white", ec=colors[s], lw=0.8, alpha=0.95), zorder=8)
draw_trajectory(ax)
ax.set_xlim(xlim); ax.set_ylim(ylim)
ax.set_xlabel("PH-only t-SNE 1"); ax.set_ylabel("PH-only t-SNE 2")
ax.set_title("PH substate structure", pad=5)
ax.set_box_aspect(1)
ax.text(0.03, 0.04, "centroid trajectory: s1 → s0 → s2 → s3", transform=ax.transAxes,
        ha="left", va="bottom", fontsize=6.2, color="#475569",
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.85))
for spine in ("top", "right"): ax.spines[spine].set_visible(False)

right = fig.add_gridspec(2, 1, left=0.55, right=0.98, bottom=0.18, top=0.88, hspace=0.32)
for j, grp in enumerate(["Control", "Hemospasia"]):
    ax = fig.add_subplot(right[j, 0])
    z = d[d["group"] == grp]
    for s in [0, 1, 2, 3]:
        zz = z[z["subcluster"] == s]
        ax.scatter(zz["tsne1"], zz["tsne2"], s=8, alpha=0.68, color=colors[s], edgecolors="none", rasterized=True)
    draw_trajectory(ax)
    ax.set_xlim(xlim); ax.set_ylim(ylim)
    ax.set_box_aspect(1)
    ax.set_title(f"{grp} (n={len(z)})", pad=3)
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ("top", "right"): ax.spines[spine].set_visible(False)

handles = [Line2D([0], [0], marker="o", linestyle="", markersize=5, markerfacecolor=colors[s], markeredgewidth=0, label=names[s]) for s in [0,1,2,3]]
fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.06, 0.02), ncol=2, frameon=False, fontsize=6.2, handletextpad=0.4, columnspacing=1.1)
fig.text(0.55, 0.91, "Group overlay on the same PH-only embedding", ha="center", va="bottom", fontsize=9, fontweight="bold")

fig.savefig(OUT / "PH_tsne_large_points_trajectory_preview.png", dpi=300, bbox_inches="tight")
plt.close(fig)
print(OUT / "PH_tsne_large_points_trajectory_preview.png")
