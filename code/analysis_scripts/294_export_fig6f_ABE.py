# -*- coding: utf-8 -*-
"""Export the two-group Hs PH versus circulating PH A/B/E composite for Fig. 6F."""

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
OUT = ROOT / "analysis_output" / "figures_pub" / "2026-10-21_Fig6F_ABE"
OUT.mkdir(parents=True, exist_ok=True)

percell = pd.read_csv(SRC / "cross_organ_percell_scores.csv.gz")
strat = pd.read_csv(SRC / "cross_organ_cycle_stratified_means.csv")
summary = pd.read_csv(SRC / "cross_organ_adjusted_summary.csv")

labels = {
    "Hs PH": "Hs PH",
    "Blood PH": "Circulating PH",
}
colors = {
    "Hs PH": "#B2182B",
    "Blood PH": "#2166AC",
}
order = ["Hs PH", "Blood PH"]
percell = percell[percell["target"].isin(order)].copy()
strat = strat[strat["target"].isin(order)].copy()
summary = summary[summary["target"].isin(order)].copy()

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

fig = plt.figure(figsize=(7.1, 5.1))
gs = fig.add_gridspec(
    2, 2,
    height_ratios=[1.05, 0.82],
    hspace=0.48,
    wspace=0.30,
    left=0.10,
    right=0.97,
    top=0.92,
    bottom=0.10,
)

# A and B: raw cell-cycle/progenitor relationship
for j, key in enumerate(order):
    ax = fig.add_subplot(gs[0, j])
    d = percell[percell["target"] == key]
    bg = d[~d["is_target"]]
    tg = d[d["is_target"]]
    ax.hexbin(
        bg["cycle_pct"], bg["progenitor_pct"],
        gridsize=34, cmap="Greys", mincnt=1, linewidths=0,
    )
    ax.scatter(
        tg["cycle_pct"], tg["progenitor_pct"],
        s=4.0, alpha=0.22, color=colors[key], edgecolors="none", rasterized=True,
    )
    sm = strat[(strat["target"] == key) & (strat["group"] == "background")]
    ax.plot(sm["mean_cycle_pct"], sm["mean_prog_pct"], color=colors[key], lw=1.2)
    row = summary[summary["target"] == key].iloc[0]
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks(np.linspace(0, 1, 6))
    ax.set_yticks(np.linspace(0, 1, 6))
    ax.set_xlabel("Cell-cycle percentile")
    if j == 0:
        ax.set_ylabel("Progenitor program percentile")
    ax.set_title(f"{labels[key]} (n={int(row['n_target']):,})", pad=5)
    ax.text(
        0.035, 0.95, "A" if j == 0 else "B",
        transform=ax.transAxes, ha="left", va="top",
        fontsize=10, fontweight="bold",
    )
    ax.tick_params(length=3, width=0.6, pad=2)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

# E: target versus background adjusted potential
ax = fig.add_subplot(gs[1, :])
positions = np.arange(len(order))
for i, key in enumerate(order):
    d = percell[percell["target"] == key]
    bg_vals = d.loc[~d["is_target"], "adjusted_progenitor_z"].dropna().values
    tg_vals = d.loc[d["is_target"], "adjusted_progenitor_z"].dropna().values
    vp_bg = ax.violinplot([bg_vals], positions=[i - 0.18], widths=0.30, showextrema=False)
    vp_tg = ax.violinplot([tg_vals], positions=[i + 0.18], widths=0.30, showextrema=False)
    for body in vp_bg["bodies"]:
        body.set_facecolor("#BDBDBD"); body.set_edgecolor("none"); body.set_alpha(0.72)
    for body in vp_tg["bodies"]:
        body.set_facecolor(colors[key]); body.set_edgecolor("none"); body.set_alpha(0.76)
    ax.plot(
        [i - 0.18, i + 0.18],
        [np.median(bg_vals), np.median(tg_vals)],
        color="#333333", lw=0.8, zorder=3,
    )
    ax.scatter(
        [i - 0.18, i + 0.18],
        [np.median(bg_vals), np.median(tg_vals)],
        s=9, color=["#777777", colors[key]], zorder=4,
    )
    row = summary[summary["target"] == key].iloc[0]
    ax.text(
        i + 0.18, 0.94,
        f"z = {row['mean_adjusted_prog_z']:.2f}",
        transform=ax.get_xaxis_transform(), ha="center", va="top",
        color=colors[key], fontsize=7,
    )
ax.axhline(0, color="#777777", lw=0.7, ls="--")
ax.set_xticks(positions)
ax.set_xticklabels([labels[k] for k in order], fontsize=8)
ax.set_ylabel("Adjusted progenitor potential (z)")
ax.set_title("Progenitor identity after cell-cycle adjustment", pad=7)
ax.text(-0.055, 1.04, "E", transform=ax.transAxes, ha="left", va="bottom",
        fontsize=10, fontweight="bold")
ax.legend(
    handles=[
        Patch(facecolor="#BDBDBD", label="background"),
        Patch(facecolor=colors["Hs PH"], label="Hs PH target"),
        Patch(facecolor=colors["Blood PH"], label="Circulating PH target"),
    ],
    frameon=False, fontsize=6.3, loc="upper left", ncol=3,
    handlelength=1.1, columnspacing=0.9,
)
ax.tick_params(length=3, width=0.6)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

stem = OUT / "Fig6F_ABE_HsPH_vs_CirculatingPH"
fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
fig.savefig(stem.with_suffix(".svg"), format="svg", bbox_inches="tight")
fig.savefig(stem.with_suffix(".pdf"), format="pdf", bbox_inches="tight")
fig.savefig(
    stem.with_suffix(".tiff"), format="tiff", dpi=1200,
    bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"},
)
plt.close(fig)

summary[summary["target"].isin(order)].to_csv(
    OUT / "Fig6F_ABE_summary.csv", index=False, encoding="utf-8-sig"
)
print(stem)
