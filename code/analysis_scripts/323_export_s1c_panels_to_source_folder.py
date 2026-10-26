# -*- coding: utf-8 -*-
"""Export the three S1-C subpanels as standalone publication figures."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager as fm
from matplotlib.patches import Rectangle

ROOT = Path(r"D:\桌面\single cell")
OUT = ROOT / "analysis_output"
FINAL = OUT / "figures_pub" / "2026-10-23_Fig1_Fig2_final" / "Supplementary_final"

HEAT_CSV = OUT / "lyo_cross_organ_mapping_scores.csv"
BEST_CSV = OUT / "Fig_lyo_cluster_best_matches.csv"
UCELL_CSV = OUT / "figures_pub" / "Supplementary_niche" / "Supp14_LyO_progenitor_Hs_UCell.csv"
OVERLAP_CSV = OUT / "figures_pub" / "Supplementary_niche" / "Supp14_LyO_progenitor_overlap.csv"

for path in (HEAT_CSV, BEST_CSV, UCELL_CSV, OVERLAP_CSV):
    if not path.exists():
        raise FileNotFoundError(path)

assert fm.findfont("Times New Roman", fallback_to_default=False)
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman"],
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 9,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "axes.linewidth": 0.65,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "legend.frameon": False,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "savefig.facecolor": "white",
})

CELL_ORDER = ["PH", "GH", "TS", "MH2", "F", "R", "B", "M", "MyoEpi", "PM-like"]


def save_panel(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), format="svg", bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), format="pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(
        stem.with_suffix(".tiff"),
        format="tiff",
        dpi=1200,
        bbox_inches="tight",
        facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)
    print(stem.name)


# S1-C1: LyO cluster-to-Hs marker mapping.
heat = pd.read_csv(HEAT_CSV)
heat = heat.loc[heat["organ"] == "Hepatic sinus"].copy()
best = pd.read_csv(BEST_CSV)
name_map = dict(zip(best["lyo_cluster"], best["lyo_name"]))
matrix = (
    heat.pivot(index="lyo_cluster", columns="target_group", values="z_score")
    .reindex(columns=CELL_ORDER)
    .sort_index()
)
labels = [f"{int(cluster)} {name_map.get(int(cluster), '')}" for cluster in matrix.index]

fig, ax = plt.subplots(figsize=(5.7, 5.6))
limit = float(np.nanmax(np.abs(matrix.values)))
image = ax.imshow(
    matrix.values,
    aspect="auto",
    cmap="RdBu_r",
    vmin=-limit,
    vmax=limit,
    interpolation="nearest",
)
ax.set_xticks(range(len(CELL_ORDER)))
ax.set_xticklabels(CELL_ORDER, rotation=55, ha="right")
ax.set_yticks(range(len(labels)))
ax.set_yticklabels(labels, fontsize=6.5)
ax.set_xlabel("Hs cell state")
ax.set_ylabel("Lymphoid-organ cluster")
ax.set_title("LyO cluster-to-Hs marker mapping", pad=5)
if 13 in matrix.index:
    row = list(matrix.index).index(13)
    ax.add_patch(Rectangle((-0.5, row - 0.5), len(CELL_ORDER), 1, fill=False, ec="#111111", lw=1))
for spine in ax.spines.values():
    spine.set_visible(False)
colorbar = fig.colorbar(image, ax=ax, fraction=0.04, pad=0.03)
colorbar.set_label("Marker-overlap z", fontsize=7)
colorbar.ax.tick_params(labelsize=6)
fig.tight_layout()
save_panel(fig, FINAL / "S1_C_1_LyO_cluster_to_Hs_marker_mapping")


# S1-C2: LyO progenitor UCell score across Hs states.
ucell = pd.read_csv(UCELL_CSV)
ucell = ucell.set_index(ucell.columns[0]).reindex(CELL_ORDER)
mean = ucell["mean"].astype(float).to_numpy()
fraction = ucell["frac_above_global"].astype(float).to_numpy()

fig, ax = plt.subplots(figsize=(5.2, 3.8))
y = np.arange(len(CELL_ORDER))
colors = ["#B2182B" if state == "PH" else "#8A8A8A" for state in CELL_ORDER]
ax.hlines(y, 0.46, mean, color="#D5D5D5", lw=1.1)
ax.scatter(mean, y, s=22 + fraction * 70, c=colors, zorder=3, edgecolors="#333333", linewidths=0.35)
for index, (value, frac) in enumerate(zip(mean, fraction)):
    ax.text(
        value + 0.003,
        index,
        f"{value:.3f} / {100 * frac:.0f}%",
        va="center",
        fontsize=6.5,
        color="#B2182B" if CELL_ORDER[index] == "PH" else "#444444",
    )
ax.set_yticks(y)
ax.set_yticklabels(CELL_ORDER)
ax.invert_yaxis()
ax.set_xlim(0.46, 0.55)
ax.set_xlabel("Mean LyO progenitor UCell score")
ax.set_title("LyO progenitor program in Hs", pad=5)
fig.tight_layout()
save_panel(fig, FINAL / "S1_C_2_LyO_progenitor_UCell")


# S1-C3: LyO progenitor top50 marker-program overlap.
overlap = pd.read_csv(OVERLAP_CSV)
overlap["label"] = overlap["label"].replace({"Niche": "Niche"})
order = ["PH", "Niche", "TS", "MH2", "MyoEpi", "M"]
overlap = overlap.set_index("label").reindex(order).reset_index()
x = np.arange(len(overlap))

fig, ax = plt.subplots(figsize=(4.6, 4.0))
ax.scatter(
    x,
    overlap["fold"],
    s=[75 if label == "PH" else 24 for label in overlap["label"]],
    c=["#B2182B" if label == "PH" else "#8A8A8A" for label in overlap["label"]],
    zorder=3,
    edgecolors="#333333",
    linewidths=0.4,
)
for index, row in overlap.iterrows():
    if row["fold"] > 0:
        ax.annotate(
            f"{int(row['overlap'])} genes; {row['fold']:.1f}×\nP={row['p_hypergeom']:.1e}",
            (index, row["fold"]),
            xytext=(8, 4),
            textcoords="offset points",
            fontsize=6.5,
            color="#B2182B",
            ha="left",
            va="bottom",
        )
    else:
        ax.text(index, 1.5, "0 genes\nP=1", ha="center", va="bottom", fontsize=6.2, color="#666666")
ax.set_xticks(x)
ax.set_xticklabels(overlap["label"], rotation=40, ha="right")
ax.set_ylim(0, 80)
ax.set_xlim(-0.6, 5.6)
ax.set_ylabel("Fold enrichment")
ax.set_title("LyO progenitor top50 overlap", pad=5)
fig.tight_layout()
save_panel(fig, FINAL / "S1_C_3_LyO_progenitor_top50_overlap")
