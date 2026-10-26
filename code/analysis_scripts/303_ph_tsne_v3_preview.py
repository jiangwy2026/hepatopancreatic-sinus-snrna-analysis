# -*- coding: utf-8 -*-
"""Render the selected smooth PH-only t-SNE parameter set."""

from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.lines import Line2D

ROOT = Path(r"D:\桌面\single cell")
IN = ROOT / "analysis_output" / "figures_pub" / "2026-10-23_Fig2D_PH_tSNE_params" / "PH_tsne_v3_p50_pc20.csv"
OUT = ROOT / "analysis_output" / "figures_pub" / "2026-10-23_Fig2D_PH_tSNE_v3"
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
    "pdf.fonttype": 42,
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

fig = plt.figure(figsize=(8.2, 4.3))
gs = fig.add_gridspec(1, 2, width_ratios=[1.35, 1.0], wspace=0.28, left=0.06, right=0.98, bottom=0.18, top=0.88)

ax = fig.add_subplot(gs[0, 0])
for s in [0, 1, 2, 3]:
    z = d[d["subcluster"] == s]
    ax.scatter(z["tsne1"], z["tsne2"], s=5.0, alpha=0.72, color=colors[s], edgecolors="none", rasterized=True)
for s in [0, 1, 2, 3]:
    z = d[d["subcluster"] == s]
    ax.text(z["tsne1"].median(), z["tsne2"].median(), f"s{s}", ha="center", va="center",
            fontsize=8, fontweight="bold", color="#1e293b",
            bbox=dict(boxstyle="round,pad=0.18", fc="white", ec=colors[s], lw=0.8, alpha=0.9))
ax.set_xlim(xlim); ax.set_ylim(ylim)
ax.set_xlabel("PH-only t-SNE 1"); ax.set_ylabel("PH-only t-SNE 2")
ax.set_title("PH substate structure", pad=5)
ax.set_box_aspect(1)
for spine in ("top", "right"): ax.spines[spine].set_visible(False)

right = fig.add_gridspec(2, 1, left=0.55, right=0.98, bottom=0.18, top=0.88, hspace=0.32)
for j, grp in enumerate(["Control", "Hemospasia"]):
    ax = fig.add_subplot(right[j, 0])
    z = d[d["group"] == grp]
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

fig.savefig(OUT / "Fig2D_PH_only_tSNE_v3_preview.png", dpi=300, bbox_inches="tight")
fig.savefig(OUT / "Fig2D_PH_only_tSNE_v3_preview.svg", format="svg", bbox_inches="tight")
fig.savefig(OUT / "Fig2D_PH_only_tSNE_v3_preview.tiff", format="tiff", dpi=1200, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
plt.close(fig)
print(OUT / "Fig2D_PH_only_tSNE_v3_preview.png")
