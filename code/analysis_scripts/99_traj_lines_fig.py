# -*- coding: utf-8 -*-
"""
99_traj_lines_fig.py - two ways of drawing lineage trajectories on the UMAP:
  A) data-driven: the minimum-spanning-tree path over cluster centroids (a single chain)
  B) hypothesis/reference style: two lineages (granulocyte: PH->GH ; monocytic: PH->MH2->MH1)
Curves are smooth splines through the cluster centroids of each path.
"""
import os, sys, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Ellipse
from matplotlib.gridspec import GridSpec
import numpy as np, pandas as pd
from scipy.interpolate import UnivariateSpline
import anndata as ad

OUT = r"D:\桌面\single cell\analysis_output"; FIG = os.path.join(OUT, "figures_pub")
sys.path.insert(0, OUT)
from pub_style import apply_pub_style
apply_pub_style(7)
ORD = ["PH", "GH", "MH1", "MH2"]
COL = {"PH": "#B2182B", "GH": "#E08214", "MH1": "#FDB863", "MH2": "#8C510A"}

def ell(pts, nsd):
    mu = pts.mean(0); cov = np.cov(pts.T)
    ev, evec = np.linalg.eigh(cov); o = ev.argsort()[::-1]; ev, evec = ev[o], evec[:, o]
    return mu, (2 * nsd * np.sqrt(ev[0]), 2 * nsd * np.sqrt(ev[1])), np.degrees(np.arctan2(evec[1, 0], evec[0, 0]))

def smooth_through(pts, n=200, s=0.0):
    t = np.arange(len(pts), dtype=float)
    fx = UnivariateSpline(t, pts[:, 0], k=min(3, len(pts) - 1), s=s)
    fy = UnivariateSpline(t, pts[:, 1], k=min(3, len(pts) - 1), s=s)
    tt = np.linspace(0, len(pts) - 1, n)
    return np.c_[fx(tt), fy(tt)]

def arrows(ax, curve, color, idx=(20, 60, 100, 140, 180), size=10):
    for i in idx:
        if i + 4 >= len(curve): continue
        ax.annotate("", xy=curve[i + 4], xytext=curve[i],
                    arrowprops=dict(arrowstyle="-|>,head_length=0.5,head_width=0.22",
                                    color=color, lw=0, shrinkA=0, shrinkB=0, mutation_scale=size), zorder=8)

def main():
    A = ad.read_h5ad(os.path.join(OUT, "adata_hemocyte_traj_geodesic.h5ad"))
    ct = A.obs["ct"].astype(str).values
    U = np.asarray(A.obsm["X_umap"]); raw = np.asarray(A.obs["geodesic_pseudotime"], dtype=float)
    root = int(json.load(open(os.path.join(OUT, "TRAJ_PSEUDOTIME_REPORT.json")))["root_cell"])
    pt = np.argsort(np.argsort(raw)) / (len(raw) - 1.0)
    cent = {g: U[ct == g].mean(0) for g in ORD}
    x0, x1 = np.percentile(U[:, 0], [0.3, 99.7]); y0, y1 = np.percentile(U[:, 1], [0.3, 99.7])
    px = 0.05 * (x1 - x0); py = 0.05 * (y1 - y0)

    variants = [("A_chain", [("PH", "MH2", "MH1", "GH")], {"PH->MH2->MH1->GH": "#7F4A12"}),
                ("B_two_lineages", [("PH", "GH"), ("PH", "MH2", "MH1")],
                 {"PH->GH": "#2E6DA4", "PH->MH2->MH1": "#7F4A12"})]
    for tag, lineages, lcol in variants:
        fig = plt.figure(figsize=(6.6, 3.5))
        gs = GridSpec(2, 2, height_ratios=[0.055, 1.0], hspace=0.06, wspace=0.10,
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
            axE.add_patch(Ellipse(mu, w, h, angle=ang, fill=False, ls=(0, (3.0, 2.4)), lw=0.55, ec="#9A9A9A", zorder=5))
        for path in lineages:
            p = np.vstack([cent[g] for g in path])
            curve = smooth_through(p)
            name = "->".join(path); c = lcol[name]
            axE.plot(curve[:, 0], curve[:, 1], lw=1.5, color=c, solid_capstyle="round", zorder=7,
                     path_effects=[pe.withStroke(linewidth=2.6, foreground="white")])
            arrows(axE, curve, c)
            i = len(curve) // 2
            axE.text(curve[i, 0], curve[i, 1] + 1.9, name, fontsize=5.6, color=c, ha="center", va="bottom",
                     zorder=9, path_effects=[pe.withStroke(linewidth=1.7, foreground="white")])
        for key, lab, col in [("GH", "Granulocyte lineage", "#2E6DA4"),
                              ("MH", "Monocytic hemocyte lineage", "#7F4A12"),
                              ("PH", "Prohemocyte", "#B2182B")]:
            m = (ct == "GH") if key == "GH" else ((ct == "MH1") | (ct == "MH2")) if key == "MH" else (ct == "PH")
            mu = U[m].mean(0)
            off = {"GH": (0.4, 4.3), "MH": (-5.4, -5.9), "PH": (1.4, 4.6)}[key]
            axE.text(mu[0] + off[0], mu[1] + off[1], lab, fontsize=5.8, color=col, ha="center", va="center",
                     zorder=9, path_effects=[pe.withStroke(linewidth=1.7, foreground="white")])
        axE.legend(loc="lower left", bbox_to_anchor=(-0.01, -0.01), ncol=2, frameon=False, markerscale=4.5,
                   handletextpad=0.25, columnspacing=0.7, fontsize=5.8, borderpad=0.1, labelspacing=0.25)
        for a in (axD, axE):
            a.set_xlim(x0 - px, x1 + px); a.set_ylim(y0 - py, y1 + py)
            a.set_xticks([]); a.set_yticks([]); a.set_xlabel("UMAP 1"); a.set_ylabel("UMAP 2")
            a.set_box_aspect(1); a.spines["top"].set_visible(True); a.spines["right"].set_visible(True)
        axD.text(-0.16, 1.06, "D.", transform=axD.transAxes, fontsize=9, fontweight="bold", va="top")
        axE.text(-0.16, 1.06, "E.", transform=axE.transAxes, fontsize=9, fontweight="bold", va="top")
        p = os.path.join(FIG, "Fig4_traj_lines_%s_preview.png" % tag)
        fig.savefig(p, dpi=400); plt.close(fig)
        print("saved", p)

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\traj_lines.log", "w", encoding="utf8", buffering=1)
    main()
