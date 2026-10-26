# -*- coding: utf-8 -*-
"""95_traj_figure_v2.py - cleaner version of the pseudotime trajectory figure (reference-style D/E)."""
import os, sys, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Ellipse
import numpy as np, pandas as pd
import anndata as ad

OUT = r"D:\桌面\single cell\analysis_output"
FIG = os.path.join(OUT, "figures_pub")
sys.path.insert(0, OUT)
from pub_style import apply_pub_style
apply_pub_style(7)
ORD = ["PH", "GH", "MH1", "MH2"]
COL = {"PH": "#B2182B", "GH": "#E08214", "MH1": "#FDB863", "MH2": "#8C510A"}

def ell(pts, nsd):
    mu = pts.mean(0); cov = np.cov(pts.T)
    ev, evec = np.linalg.eigh(cov); o = ev.argsort()[::-1]
    ev, evec = ev[o], evec[:, o]
    ang = np.degrees(np.arctan2(evec[1, 0], evec[0, 0]))
    return mu, (2 * nsd * np.sqrt(ev[0]), 2 * nsd * np.sqrt(ev[1])), ang

def main():
    A = ad.read_h5ad(os.path.join(OUT, "adata_hemocyte_traj_geodesic.h5ad"))
    ct = A.obs["ct"].astype(str).values
    U = np.asarray(A.obsm["X_umap"]); pt = np.asarray(A.obs["geodesic_pseudotime"], dtype=float)
    root = int(json.load(open(os.path.join(OUT, "TRAJ_PSEUDOTIME_REPORT.json")))["root_cell"])
    x0, x1 = np.percentile(U[:, 0], [0.4, 99.6]); y0, y1 = np.percentile(U[:, 1], [0.4, 99.6])
    padx = 0.05 * (x1 - x0); pady = 0.05 * (y1 - y0)

    fig, axes = plt.subplots(1, 2, figsize=(6.6, 3.2))
    # ---- D: pseudotime
    ax = axes[0]
    o = np.argsort(-pt)
    s = ax.scatter(U[o, 0], U[o, 1], c=pt[o], s=1.0, cmap="Blues_r", linewidths=0,
                   rasterized=True, vmin=0, vmax=1)
    ax.scatter(U[root, 0], U[root, 1], s=22, marker="o", facecolor="none", edgecolor="#B2182B", lw=0.8, zorder=6)
    ax.text(U[root, 0] - 1.2, U[root, 1] + 3.4, "root (PH)", fontsize=6, color="#B2182B", ha="center",
            path_effects=[pe.withStroke(linewidth=1.6, foreground="white")], zorder=7)
    cb = fig.colorbar(s, ax=ax, orientation="horizontal", location="top", pad=0.015, fraction=0.045, aspect=30)
    cb.set_label("Pseudotime", fontsize=6.5); cb.ax.tick_params(labelsize=5.5, length=1.4, pad=1)
    cb.outline.set_visible(False)
    ax.set_title("Prohemocyte-rooted pseudotime", fontsize=7, pad=12)
    # ---- E: groups + lineages
    ax = axes[1]
    for g in ORD:
        m = ct == g
        ax.scatter(U[m, 0], U[m, 1], s=1.0, c=COL[g], linewidths=0, rasterized=True, label="%s (%d)" % (g, int(m.sum())))
    spec = [("GH", (ct == "GH"), "#2E6DA4", "Granulocyte lineage", (-0.5, 4.6)),
            ("MH", (ct == "MH1") | (ct == "MH2"), "#7F4A12", "Monocytic hemocyte lineage", (-4.6, -5.2)),
            ("PH", (ct == "PH"), "#B2182B", "Prohemocyte", (2.6, -3.6))]
    for _, m, c, lab, off in spec:
        mu, (w, h), ang = ell(U[m], 1.7)
        ax.add_patch(Ellipse(mu, w, h, angle=ang, fill=False, ls=(0, (3.4, 2.0)), lw=0.7, ec=c, zorder=5))
        ax.text(mu[0] + off[0], mu[1] + off[1], lab, fontsize=5.8, color=c, ha="center", va="center", zorder=7,
                path_effects=[pe.withStroke(linewidth=1.7, foreground="white")])
    ax.legend(loc="lower left", bbox_to_anchor=(0.0, 0.0), ncol=2, frameon=False, markerscale=4.5,
              handletextpad=0.25, columnspacing=0.7, fontsize=5.8, borderpad=0.1, labelspacing=0.25)
    ax.set_title("Functional groups and lineages", fontsize=7, pad=12)
    for a in axes:
        a.set_xlim(x0 - padx, x1 + padx); a.set_ylim(y0 - pady, y1 + pady)
        a.set_xticks([]); a.set_yticks([])
        a.set_xlabel("UMAP 1"); a.set_ylabel("UMAP 2")
    for a, L in zip(axes, ["D.", "E."]):
        a.text(-0.10, 1.13, L, transform=a.transAxes, fontsize=9, fontweight="bold", va="top")
    fig.subplots_adjust(wspace=0.14)
    p = os.path.join(FIG, "Fig4_traj_pseudotime_preview.png")
    fig.savefig(p, dpi=400, bbox_inches="tight")
    print("saved", p, flush=True)

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\traj_figure.log", "w", encoding="utf8", buffering=1)
    main()
