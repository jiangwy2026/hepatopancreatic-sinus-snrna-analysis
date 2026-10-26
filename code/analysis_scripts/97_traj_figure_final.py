# -*- coding: utf-8 -*-
"""97_traj_figure_final.py - reference-style trajectory figure (final preview, Fig4 D/E)."""
import os, sys, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Ellipse
from matplotlib.gridspec import GridSpec
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
    U = np.asarray(A.obsm["X_umap"]); raw = np.asarray(A.obs["geodesic_pseudotime"], dtype=float)
    root = int(json.load(open(os.path.join(OUT, "TRAJ_PSEUDOTIME_REPORT.json")))["root_cell"])
    pt = np.argsort(np.argsort(raw)) / (len(raw) - 1.0)
    x0, x1 = np.percentile(U[:, 0], [0.3, 99.7]); y0, y1 = np.percentile(U[:, 1], [0.3, 99.7])
    px = 0.05 * (x1 - x0); py = 0.05 * (y1 - y0)

    fig = plt.figure(figsize=(6.6, 3.5))
    gs = GridSpec(2, 2, height_ratios=[0.055, 1.0], width_ratios=[1, 1], hspace=0.06, wspace=0.10,
                  left=0.055, right=0.985, top=0.90, bottom=0.09)
    axc = fig.add_subplot(gs[0, 0]); axD = fig.add_subplot(gs[1, 0]); axE = fig.add_subplot(gs[1, 1])

    o = np.argsort(-pt)
    s = axD.scatter(U[o, 0], U[o, 1], c=pt[o], s=1.0, cmap="Blues_r", linewidths=0, rasterized=True, vmin=0, vmax=1)
    axD.scatter(U[root, 0], U[root, 1], s=20, facecolor="none", edgecolor="#B2182B", lw=0.8, zorder=6)
    axD.text(U[root, 0] - 1.4, U[root, 1] + 3.2, "root (PH)", fontsize=6, color="#B2182B", ha="center",
             path_effects=[pe.withStroke(linewidth=1.6, foreground="white")], zorder=7)
    cb = fig.colorbar(s, cax=axc, orientation="horizontal")
    cb.ax.xaxis.set_ticks_position("top"); cb.ax.xaxis.set_label_position("top")
    cb.set_label("Pseudotime", fontsize=6.5, labelpad=1.5)
    cb.ax.tick_params(labelsize=5.5, length=1.4, pad=1); cb.outline.set_visible(False)

    for g in ORD:
        m = ct == g
        axE.scatter(U[m, 0], U[m, 1], s=1.0, c=COL[g], linewidths=0, rasterized=True,
                    label="%s (%d)" % (g, int(m.sum())))
    for m, c, lab, off in [((ct == "GH"), "#2E6DA4", "Granulocyte lineage", (0.4, 4.3)),
                           ((ct == "MH1") | (ct == "MH2"), "#7F4A12", "Monocytic hemocyte lineage", (-5.4, -5.9)),
                           ((ct == "PH"), "#B2182B", "Prohemocyte", (1.4, 4.6))]:
        mu, (w, h), ang = ell(U[m], 1.7)
        axE.add_patch(Ellipse(mu, w, h, angle=ang, fill=False, ls=(0, (3.4, 2.0)), lw=0.7, ec=c, zorder=5))
        axE.text(mu[0] + off[0], mu[1] + off[1], lab, fontsize=5.8, color=c, ha="center", va="center", zorder=7,
                 path_effects=[pe.withStroke(linewidth=1.7, foreground="white")])
    axE.legend(loc="lower left", bbox_to_anchor=(-0.01, -0.01), ncol=2, frameon=False, markerscale=4.5,
               handletextpad=0.25, columnspacing=0.7, fontsize=5.8, borderpad=0.1, labelspacing=0.25)

    for a in (axD, axE):
        a.set_xlim(x0 - px, x1 + px); a.set_ylim(y0 - py, y1 + py)
        a.set_xticks([]); a.set_yticks([])
        a.set_xlabel("UMAP 1"); a.set_ylabel("UMAP 2")
        a.set_box_aspect(1)
        a.spines["top"].set_visible(True); a.spines["right"].set_visible(True)
    axD.text(-0.16, 1.06, "D.", transform=axD.transAxes, fontsize=9, fontweight="bold", va="top")
    axE.text(-0.16, 1.06, "E.", transform=axE.transAxes, fontsize=9, fontweight="bold", va="top")
    p = os.path.join(FIG, "Fig4_traj_pseudotime_preview.png")
    fig.savefig(p, dpi=400)
    print("saved", p)

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\traj_figure.log", "w", encoding="utf8", buffering=1)
    main()
