# -*- coding: utf-8 -*-
"""Export a square, standalone cell-cycle-adjusted progenitor violin plot."""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import Patch

ROOT = Path(r"D:\桌面\single cell")
SRC = ROOT / "analysis_output" / "cross_organ_progenitor"
OUT = ROOT / "analysis_output" / "figures_pub" / "2026-10-22_Fig6F_adjusted_violin_square"
OUT.mkdir(parents=True, exist_ok=True)

percell = pd.read_csv(SRC / "cross_organ_percell_scores.csv.gz")
summary = pd.read_csv(SRC / "cross_organ_adjusted_summary.csv")
order = ["Hs PH", "Blood PH"]
labels = {"Hs PH": "Hs PH", "Blood PH": "Circulating PH"}
colors = {"Hs PH": "#B2182B", "Blood PH": "#2166AC"}
percell = percell[percell["target"].isin(order)].copy()
summary = summary[summary["target"].isin(order)].copy()

assert fm.findfont("Times New Roman", fallback_to_default=False)
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman"],
    "font.size": 8.5,
    "axes.titlesize": 9.5,
    "axes.labelsize": 8.5,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.linewidth": 0.7,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "savefig.facecolor": "white",
})

fig, ax = plt.subplots(figsize=(5.4, 5.4))
fig.subplots_adjust(left=0.17, right=0.96, bottom=0.13, top=0.84)

positions = np.arange(len(order))
offset = 0.23
for i, key in enumerate(order):
    d = percell[percell["target"] == key]
    bg_vals = d.loc[~d["is_target"], "adjusted_progenitor_z"].dropna().values
    tg_vals = d.loc[d["is_target"], "adjusted_progenitor_z"].dropna().values

    vp_bg = ax.violinplot(
        [bg_vals], positions=[i - offset], widths=0.50,
        bw_method=0.35, points=220, showextrema=False,
    )
    vp_tg = ax.violinplot(
        [tg_vals], positions=[i + offset], widths=0.50,
        bw_method=0.35, points=220, showextrema=False,
    )
    for body in vp_bg["bodies"]:
        body.set_facecolor("#BDBDBD")
        body.set_edgecolor("none")
        body.set_alpha(0.72)
    for body in vp_tg["bodies"]:
        body.set_facecolor(colors[key])
        body.set_edgecolor("none")
        body.set_alpha(0.76)

    ax.plot(
        [i - offset, i + offset],
        [np.median(bg_vals), np.median(tg_vals)],
        color="#333333", lw=0.9, zorder=3,
    )
    ax.scatter(
        [i - offset, i + offset],
        [np.median(bg_vals), np.median(tg_vals)],
        s=11, color=["#777777", colors[key]], zorder=4,
    )

    row = summary[summary["target"] == key].iloc[0]
    ax.text(
        i + offset, 0.82,
        f"z = {row['mean_adjusted_prog_z']:.2f}",
        transform=ax.get_xaxis_transform(),
        ha="center", va="top", color=colors[key], fontsize=7.5,
    )

ax.axhline(0, color="#777777", lw=0.8, ls="--")
ax.set_ylim(-2.6, 2.6)
ax.set_xticks(positions)
ax.set_xticklabels([labels[k] for k in order], fontsize=8.5)
ax.set_ylabel("Adjusted progenitor potential (z)")
ax.set_title("Progenitor identity after cell-cycle adjustment", pad=28)
ax.tick_params(length=3, width=0.7)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

ax.legend(
    handles=[
        Patch(facecolor="#BDBDBD", label="background"),
        Patch(facecolor=colors["Hs PH"], label="Hs PH target"),
        Patch(facecolor=colors["Blood PH"], label="Circulating PH target"),
    ],
    frameon=False,
    fontsize=6.8,
    loc="lower center",
    bbox_to_anchor=(0.5, 1.01),
    ncol=3,
    handlelength=1.1,
    columnspacing=1.0,
)

stem = OUT / "Adjusted_progenitor_identity_square"
fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches=None)
fig.savefig(stem.with_suffix(".svg"), format="svg", bbox_inches=None)
fig.savefig(stem.with_suffix(".pdf"), format="pdf", bbox_inches=None)
fig.savefig(
    stem.with_suffix(".tiff"), format="tiff", dpi=1200,
    bbox_inches=None, pil_kwargs={"compression": "tiff_lzw"},
)
plt.close(fig)
summary.to_csv(OUT / "Adjusted_violin_square_summary.csv", index=False, encoding="utf-8-sig")
print(stem)
