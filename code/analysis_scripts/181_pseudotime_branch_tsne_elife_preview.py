# -*- coding: utf-8 -*-
"""181_pseudotime_branch_tsne_elife_preview.py
eLife-style branched pseudotime preview on the all-cell t-SNE:
PH -> TS -> GH and PH -> TS -> MH2.
All ten subtypes are retained; non-hemocyte populations are grey.
"""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import FancyArrowPatch
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
C_TRUNK = "#222222"; C_GRAN = "#2166AC"; C_MONO = "#8C510A"


def normrank(x):
    x = np.asarray(x, float); good = np.isfinite(x); out = np.full(len(x), np.nan)
    r = rankdata(x[good], method="average")
    out[good] = (r - 1) / max(1, good.sum() - 1)
    return out


def main():
    T = ad.read_h5ad(os.path.join(OUT, "adata_soupx_hvg_tsne.h5ad"))
    Hh = ad.read_h5ad(os.path.join(OUT, "adata_hemocyte_traj_geodesic.h5ad"))
    U = np.asarray(T.obsm["X_tsne"])
    ct = T.obs["vendor_cluster"].astype(str).map(GROUPS).values
    dpt = np.asarray(Hh.obs["dpt_pseudotime"], float); dpt[~np.isfinite(dpt)] = np.nan
    pt = pd.Series(normrank(dpt), index=Hh.obs_names).reindex(T.obs_names).values
    hem = np.isfinite(pt)
    fig = plt.figure(figsize=(4.45, 4.05), facecolor="white")
    ax = fig.add_axes([0.10, 0.10, 0.73, 0.78])
    ax.scatter(U[~hem, 0], U[~hem, 1], s=1.0, c="#D7D7D7", linewidths=0,
               rasterized=True, zorder=1)
    idx = np.where(hem)[0][np.argsort(-pt[hem])]
    s = ax.scatter(U[idx, 0], U[idx, 1], c=pt[idx], s=1.1, cmap="viridis",
                   linewidths=0, rasterized=True, vmin=0, vmax=1, zorder=2)
    cent = {g: U[ct == g].mean(0) for g in ORDER if (ct == g).sum() > 0}
    # Shared trunk and two branches.
    branches = [
        (cent["PH"], cent["TS"], 0.02, C_TRUNK, 2.2),
        (cent["TS"], cent["GH"], 0.25, C_GRAN, 1.8),
        (cent["TS"], cent["MH2"], -0.22, C_MONO, 1.8),
    ]
    for p1, p2, rad, col, lw in branches:
        ax.add_patch(FancyArrowPatch(p1, p2, connectionstyle="arc3,rad=%s" % rad,
                                     arrowstyle="-|>,head_length=0.62,head_width=0.30",
                                     mutation_scale=12, lw=lw, color=col, alpha=0.98,
                                     zorder=12, shrinkA=3, shrinkB=3,
                                     path_effects=[pe.withStroke(linewidth=3.2, foreground="white")]))
    # Branch labels
    def mid(a, b, frac=0.55):
        p = a + frac * (b - a)
        return p[0], p[1]
    x, y = mid(cent["PH"], cent["TS"], 0.55)
    ax.text(x - 0.8, y + 0.45, "PH → TS", fontsize=5.2, color=C_TRUNK, ha="center", va="center",
            path_effects=[pe.withStroke(linewidth=1.8, foreground="white")], zorder=14)
    x, y = mid(cent["TS"], cent["GH"], 0.62)
    ax.text(x, y - 0.5, "granular branch", fontsize=5.1, color=C_GRAN, ha="center", va="top",
            path_effects=[pe.withStroke(linewidth=1.8, foreground="white")], zorder=14)
    x, y = mid(cent["TS"], cent["MH2"], 0.62)
    ax.text(x, y + 0.4, "phagocytic branch", fontsize=5.1, color=C_MONO, ha="center", va="bottom",
            path_effects=[pe.withStroke(linewidth=1.8, foreground="white")], zorder=14)
    # labels for all ten states
    offs = {"PH": (-18, 8), "GH": (-20, -11), "TS": (0, -15), "MH2": (0, 8),
            "F": (0, -9), "R": (12, 0), "B": (8, 4), "M": (-12, -4),
            "MyoEpi": (10, -7), "PM-like": (0, 6)}
    for g in ORDER:
        if g not in cent: continue
        c = cent[g]; is_hem = g in HEM
        ax.annotate(g, xy=(c[0], c[1]), xytext=offs[g], textcoords="offset points",
                    fontsize=5.7, weight="bold", ha="center", va="center",
                    color="white" if is_hem else "#444444",
                    path_effects=[pe.withStroke(linewidth=1.9, foreground="#333333" if is_hem else "white")],
                    zorder=15)
    cax = fig.add_axes([0.855, 0.56, 0.020, 0.23])
    cb = fig.colorbar(s, cax=cax)
    cb.set_label("DPT", fontsize=5.0, labelpad=1.5)
    cb.ax.tick_params(labelsize=5.0, length=0.9, pad=1.0)
    cb.set_ticks([0, 0.5, 1.0]); cb.outline.set_visible(False)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("t-SNE 1", fontsize=6.0); ax.set_ylabel("t-SNE 2", fontsize=6.0)
    ax.set_box_aspect(1)
    for sp in ax.spines.values():
        sp.set_color("#777777"); sp.set_linewidth(0.45)
    ax.set_title("Branched pseudotime: PH → TS → GH / MH2", fontsize=6.8, pad=5)
    fig.text(0.50, 0.018, "Grey: non-hemocyte populations. Branches are DPT-based candidate paths.",
             ha="center", va="bottom", fontsize=5.0, color="#333333")
    p = os.path.join(PREV, "02c_pseudotime_branch_tsne_preview.png")
    fig.savefig(p, dpi=400, bbox_inches="tight", pad_inches=0.04, facecolor="white")
    print("saved", p)
    plt.close(fig)

if __name__ == "__main__":
    main()
