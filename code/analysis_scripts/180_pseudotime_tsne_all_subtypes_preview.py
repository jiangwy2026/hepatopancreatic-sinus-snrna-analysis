# -*- coding: utf-8 -*-
"""180_pseudotime_tsne_all_subtypes_preview.py
PH-rooted DPT shown on the all-cell t-SNE. Hemocyte states retain the
pseudotime gradient; all other populations remain visible in grey and labelled.
Standalone preview.
"""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import FancyArrowPatch
from matplotlib.colors import Normalize
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
from scipy.stats import rankdata

sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"
PREV = os.path.join(OUT, "figures_pub", "sc_summary_new_names", "individual_previews")
os.makedirs(PREV, exist_ok=True)
GROUPS = {"0": "TS", "1": "GH", "2": "PH", "3": "MH2", "4": "PM-like",
          "5": "F", "6": "R", "7": "B", "8": "MyoEpi", "9": "M"}
HEM = ["PH", "GH", "TS", "MH2"]
ORDER = ["PH", "GH", "TS", "MH2", "F", "R", "B", "M", "MyoEpi", "PM-like"]


def normrank(x):
    x = np.asarray(x, float); good = np.isfinite(x); out = np.full(len(x), np.nan)
    if good.sum():
        r = rankdata(x[good], method="average")
        out[good] = (r - 1) / max(1, good.sum() - 1)
    return out


def main():
    T = ad.read_h5ad(os.path.join(OUT, "adata_soupx_hvg_tsne.h5ad"))
    Hh = ad.read_h5ad(os.path.join(OUT, "adata_hemocyte_traj_geodesic.h5ad"))
    U = np.asarray(T.obsm["X_tsne"])
    ct = T.obs["vendor_cluster"].astype(str).map(GROUPS).values
    dpt = np.asarray(Hh.obs["dpt_pseudotime"], float)
    dpt[~np.isfinite(dpt)] = np.nan
    ptr = pd.Series(normrank(dpt), index=Hh.obs_names)
    pt_all = ptr.reindex(T.obs_names).values
    hem_mask = np.isfinite(pt_all)

    fig = plt.figure(figsize=(4.15, 3.90), facecolor="white")
    ax = fig.add_axes([0.115, 0.105, 0.74, 0.78])
    # non-hemocyte populations: grey but retained
    other = ~hem_mask
    ax.scatter(U[other, 0], U[other, 1], s=1.05, c="#D5D5D5", linewidths=0,
               rasterized=True, zorder=1)
    # hemocyte populations: PH-rooted DPT gradient
    o = np.argsort(-pt_all[hem_mask])
    idx = np.where(hem_mask)[0][o]
    sc0 = ax.scatter(U[idx, 0], U[idx, 1], c=pt_all[idx], s=1.15, cmap="viridis",
                     linewidths=0, rasterized=True, vmin=0, vmax=1, zorder=2)
    cent = {g: U[ct == g].mean(0) for g in ORDER if (ct == g).sum()}
    # candidate paths based on DPT ordering
    paths = [(cent["PH"], cent["TS"], 0.04, "#111111", 1.9),
             (cent["TS"], cent["GH"], 0.12, "#2166AC", 1.5),
             (cent["TS"], cent["MH2"], -0.12, "#8C510A", 1.5)]
    for p1, p2, rad, col, lw in paths:
        ax.add_patch(FancyArrowPatch(p1, p2, connectionstyle="arc3,rad=%s" % rad,
                                     arrowstyle="-|>,head_length=0.62,head_width=0.30",
                                     mutation_scale=11, lw=lw, color=col, alpha=0.98,
                                     zorder=12, shrinkA=3, shrinkB=3,
                                     path_effects=[pe.withStroke(linewidth=3.0, foreground="white")]))
    # PH root is indicated by the PH label; no single-cell star is shown.
    # labels: all ten populations
    offs = {"PH": (-20, 9), "GH": (-18, -11), "TS": (0, -15), "MH2": (0, 8),
            "F": (0, -9), "R": (12, 0), "B": (8, 4), "M": (-12, -4),
            "MyoEpi": (10, -7), "PM-like": (0, 6)}
    for g in ORDER:
        if g not in cent: continue
        c = cent[g]
        is_hem = g in HEM
        lab = g
        ax.annotate(lab, xy=(c[0], c[1]), xytext=offs[g], textcoords="offset points",
                    fontsize=5.7, weight="bold", ha="center", va="center",
                    color="white" if is_hem else "#444444",
                    path_effects=[pe.withStroke(linewidth=1.9, foreground="#333333" if is_hem else "white")],
                    zorder=10)
    # smaller inset colour bar
    cax = fig.add_axes([0.875, 0.54, 0.020, 0.23])
    cb = fig.colorbar(sc0, cax=cax)
    cb.set_label("DPT", fontsize=5.0, labelpad=1.5)
    cb.ax.tick_params(labelsize=5.0, length=0.9, pad=1.0)
    cb.set_ticks([0, 0.5, 1.0]); cb.outline.set_visible(False)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("t-SNE 1", fontsize=6.0); ax.set_ylabel("t-SNE 2", fontsize=6.0)
    ax.set_box_aspect(1)
    for sp in ax.spines.values():
        sp.set_color("#777777"); sp.set_linewidth(0.45)
    ax.set_title("Pseudotime (all subtypes)", fontsize=6.9, pad=5)
    fig.text(0.50, 0.018, "Grey: non-hemocyte populations. Arrows: DPT-based candidate paths.",
             ha="center", va="bottom", fontsize=5.0, color="#333333")
    p = os.path.join(PREV, "02b_pseudotime_tsne_all_subtypes_preview.png")
    fig.savefig(p, dpi=400, bbox_inches="tight", pad_inches=0.04, facecolor="white")
    print("saved", p)
    plt.close(fig)

if __name__ == "__main__":
    main()
