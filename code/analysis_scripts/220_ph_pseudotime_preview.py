# -*- coding: utf-8 -*-
"""220_ph_pseudotime_preview.py
PH-specific pseudotime previews:
1) ridge-density plot of PH substates by group;
2) conventional UMAP pseudotime plot with PH substate labels and velocity arrows.
"""
import os
os.environ.setdefault("MPLCONFIGDIR", r"D:\桌面\single cell\analysis_output\.mplconfig")
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

ROOT = r"D:\桌面\single cell"
OUT = os.path.join(ROOT, "analysis_output")
PARTS = os.path.join(OUT, "figures_pub", "fig5_parts")
os.makedirs(PARTS, exist_ok=True)
sys.path.insert(0, OUT)
from pub_style import apply_pub_style

apply_pub_style(7)
mpl.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Nimbus Roman", "Liberation Serif", "DejaVu Serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})
sc.settings.verbosity = 0

PH_KEYS = ["PH_s0", "PH_s1", "PH_s2", "PH_s3"]
PH_LABELS = {"PH_s0": "s0 ECM-high", "PH_s1": "s1 progenitor",
             "PH_s2": "s2 immune-primed", "PH_s3": "s3 TGM1+"}
COL = {"PH_s0": "#59A14F", "PH_s1": "#E15759", "PH_s2": "#F28E2B", "PH_s3": "#4E79A7"}


def main():
    A = sc.read_h5ad(os.path.join(OUT, "adata_velocity_processed.h5ad"))
    A.obs["parent"] = A.obs["vendor_cluster"].astype(str).map(
        {"0": "TS", "1": "GH", "2": "PH", "3": "MH2", "4": "PM-like",
         "5": "F", "6": "R", "7": "B", "8": "MyoEpi", "9": "M"})
    A.obs["group2"] = A.obs["sample"].astype(str).map({"blank": "Control", "treated": "Hemospasia"})
    sub = pd.read_csv(os.path.join(OUT, "dynamic_subpop_subclusters", "PH_subcluster_umap.csv"), index_col=0)
    sub["sub"] = "PH_s" + sub["subcluster"].astype(str)
    A.obs["sub"] = pd.Series(index=A.obs_names, dtype=object)
    A.obs.loc[sub.index, "sub"] = "PH_s" + sub["subcluster"].astype(str)
    ph = A[A.obs["parent"] == "PH"].copy()
    print("PH cells", ph.n_obs, ph.obs["sub"].value_counts().to_dict())

    # 1) ridge-density plot
    fig = plt.figure(figsize=(4.8, 3.5), facecolor="white")
    ax = fig.add_axes([0.27, 0.17, 0.68, 0.70])
    xmin, xmax = ph.obs["velocity_pseudotime"].quantile([0.01, 0.99]).astype(float).values
    pad = 0.12 * (xmax - xmin)
    xx = np.linspace(xmin - pad, xmax + pad, 240)
    summ = pd.read_csv(os.path.join(OUT, "Fig5_dynamic_subcluster_pseudotime_summary.csv"))
    sig = set(summ.loc[(summ.key.isin(PH_KEYS)) & (summ.fdr < 0.05), "key"])
    for row, key in enumerate(PH_KEYS[::-1]):
        d = ph[ph.obs["sub"].astype(str) == key]
        for g, ls, alpha in [("Control", (0, (3, 2)), 0.80), ("Hemospasia", "-", 1.0)]:
            x = d.obs.loc[d.obs.group2 == g, "velocity_pseudotime"].dropna().astype(float).values
            if len(x) >= 5:
                y = gaussian_kde(x, bw_method="scott")(xx)
                y = y / max(y.max(), 1e-9) * 0.80
                ax.plot(xx, row + y, color=COL[key], lw=0.95, ls=ls, alpha=alpha)
                ax.fill_between(xx, row, row + y, color=COL[key], alpha=0.08 if g == "Hemospasia" else 0)
                ax.plot([np.median(x), np.median(x)], [row - 0.05, row + 0.10],
                        color=COL[key], lw=1.0 if g == "Hemospasia" else 0.7, ls=ls)
    ax.set_yticks(np.arange(4))
    ax.set_yticklabels([PH_LABELS[k] + (" *" if k in sig else "") for k in PH_KEYS[::-1]], fontsize=5.2)
    ax.set_xlim(xx.min(), xx.max()); ax.set_ylim(-0.35, 3.35)
    ax.set_xlabel("Velocity pseudotime", fontsize=5.8)
    ax.tick_params(axis="x", labelsize=5.1, length=1.4, width=0.4)
    ax.tick_params(axis="y", length=0, pad=2)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.grid(axis="x", color="#EFEFEF", lw=0.35, zorder=0)
    ax.set_title("PH substate pseudotime", fontsize=7.0, pad=5)
    fig.text(0.5, 0.02, "Solid = Hemospasia; dashed = Control; vertical ticks = medians. "
             "* FDR < 0.05.", ha="center", va="bottom", fontsize=5.0, color="#555555")
    out1 = os.path.join(PARTS, "Fig5_PH_pseudotime_density_preview.png")
    fig.savefig(out1, dpi=450, bbox_inches="tight", pad_inches=0.04, facecolor="white")
    print("saved", out1)
    plt.close(fig)

    # 2) conventional global UMAP coloured by pseudotime
    U = np.asarray(ph.obsm["X_umap"])
    pt = ph.obs["velocity_pseudotime"].astype(float).values
    cl = ph.obs["sub"].astype(str).values
    fig = plt.figure(figsize=(4.6, 4.0), facecolor="white")
    ax = fig.add_axes([0.07, 0.10, 0.76, 0.80])
    sca = ax.scatter(U[:, 0], U[:, 1], c=pt, s=3.2, cmap="viridis",
                     linewidths=0, rasterized=True, alpha=0.90)
    for key in PH_KEYS:
        m = cl == key
        if m.sum() >= 3:
            ax.text(np.median(U[m, 0]), np.median(U[m, 1]), key.split("_")[1],
                    fontsize=5.2, ha="center", va="center", color="white", fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.12", fc="#333333", ec="none", alpha=0.72),
                    zorder=8)
    ir = int(np.nanargmin(pt))
    ax.scatter([U[ir, 0]], [U[ir, 1]], marker="*", s=95,
               c="#FFD92F", edgecolors="#222222", linewidths=0.5, zorder=10)
    ax.annotate("root", xy=(U[ir, 0], U[ir, 1]), xytext=(5, -9),
                textcoords="offset points", fontsize=5.0, color="#222222",
                ha="left", va="top")
    x1, x2 = np.quantile(U[:, 0], [0.01, 0.99])
    y1, y2 = np.quantile(U[:, 1], [0.01, 0.99])
    xp = 0.08 * (x2 - x1); yp = 0.08 * (y2 - y1)
    ax.set_xlim(x1 - xp, x2 + xp); ax.set_ylim(y1 - yp, y2 + yp)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_linewidth(0.4); sp.set_color("#777777")
    ax.set_title("PH pseudotime on global UMAP", fontsize=7.0, pad=5)
    cax = fig.add_axes([0.86, 0.20, 0.020, 0.58])
    cb = fig.colorbar(sca, cax=cax); cb.set_label("Velocity pseudotime", fontsize=5.5, labelpad=4)
    cb.ax.tick_params(labelsize=5.0, length=1.3, width=0.4); cb.outline.set_linewidth(0.4)
    fig.text(0.5, 0.02, "Global UMAP from the RNA velocity object; axis trimmed to 1st–99th percentile. "
             "Root = earliest PH cell.", ha="center", va="bottom", fontsize=5.0, color="#555555")
    out2 = os.path.join(PARTS, "Fig5_PH_pseudotime_umap_preview.png")
    fig.savefig(out2, dpi=450, bbox_inches="tight", pad_inches=0.04, facecolor="white")
    print("saved", out2)
    plt.close(fig)

if __name__ == "__main__":
    main()
