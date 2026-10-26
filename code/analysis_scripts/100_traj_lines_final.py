# -*- coding: utf-8 -*-
"""100_traj_lines_final.py - final preview: pseudotime (D) + groups, lineage ellipses and
trajectory curves (E).  Curves = smooth splines through the cluster centroids of each lineage,
direction given by the PH root and pseudotime."""
import os, sys, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Ellipse
from matplotlib.gridspec import GridSpec
import numpy as np
from scipy.interpolate import UnivariateSpline
import anndata as ad

OUT = r"D:\桌面\single cell\analysis_output"; FIG = os.path.join(OUT, "figures_pub")
sys.path.insert(0, OUT)
from pub_style import apply_pub_style
apply_pub_style(7)
ORD = ["PH", "GH", "MH1", "MH2"]
COL = {"PH": "#B2182B", "GH": "#E08214", "MH1": "#FDB863", "MH2": "#8C510A"}
LIN = [("PH", "GH", "#2E6DA4"), ("PH", "MH2", "MH1", "#7F4A12")]

def ell(pts, nsd):
    mu = pts.mean(0); cov = np.cov(pts.T)
    ev, evec = np.linalg.eigh(cov); o = ev.argsort()[::-1]; ev, evec = ev[o], evec[:, o]
    return mu, (2 * nsd * np.sqrt(ev[0]), 2 * nsd * np.sqrt(ev[1])), np.degrees(np.arctan2(evec[1, 0], evec[0, 0]))

def main():
    A = ad.read_h5ad(os.path.join(OUT, "adata_hemocyte_traj_geodesic.h5ad"))
    ct = A.obs["ct"].astype(str).values
    U = np.asarray(A.obsm["X_umap"]); raw = np.asarray(A.obs["geodesic_pseudotime"], dtype=float)
    root = int(json.load(open(os.path.join(OUT, "TRAJ_PSEUDOTIME_REPORT.json")))["root_cell"])
    pt = np.argsort(np.argsort(raw)) / (len(raw) - 1.0)
    cent = {g: U[ct == g].mean(0) for g in ORD}
    x0, x1 = np.percentile(U[:, 0], [0.3, 99.7]); y0, y1 = np.percentile(U[:, 1], [0.3, 99.7])
    px = 0.05 * (x1 - x0); py = 0.05 * (y1 - y0)

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
    for m, c, off in [((ct == "GH"), "#2E6DA4", (0.4, 4.3)),
                      ((ct == "MH1") | (ct == "MH2"), "#7F4A12", (-5.4, -5.9)),
                      ((ct == "PH"), "#B2182B", (1.4, 4.6))]:
        mu, (w, h), ang = ell(U[m], 1.7)
        axE.add_patch(Ellipse(mu, w, h, angle=ang, fill=False, ls=(0, (3.0, 2.4)), lw=0.55, ec="#9A9A9A", zorder=5))
    for spec in LIN:
        path = spec[:-1]; c = spec[-1]
        p = np.vstack([cent[g] for g in path])
        t = np.arange(len(p), dtype=float)
        fx = UnivariateSpline(t, p[:, 0], k=min(3, len(p) - 1)); fy = UnivariateSpline(t, p[:, 1], k=min(3, len(p) - 1))
        tt = np.linspace(0, len(p) - 1, 220); curve = np.c_[fx(tt), fy(tt)]
        axE.plot(curve[:, 0], curve[:, 1], lw=1.5, color=c, solid_capstyle="round", zorder=8,
                 path_effects=[pe.withStroke(linewidth=2.7, foreground="white")])
        for i in (len(curve) // 2, len(curve) - 6):
            if i + 4 >= len(curve): i = len(curve) - 8
            axE.annotate("", xy=curve[i + 4], xytext=curve[i],
                         arrowprops=dict(arrowstyle="-|>,head_length=0.55,head_width=0.24", color=c, lw=0,
                                         shrinkA=0, shrinkB=0, mutation_scale=11), zorder=9)
    for key, lab, col, off in [("GH", "Granulocyte lineage", "#2E6DA4", (0.4, 4.3)),
                               ("MH", "Monocytic hemocyte lineage", "#7F4A12", (-5.4, -5.9)),
                               ("PH", "Prohemocyte", "#B2182B", (1.4, 4.6))]:
        m = (ct == "GH") if key == "GH" else ((ct == "MH1") | (ct == "MH2")) if key == "MH" else (ct == "PH")
        mu = U[m].mean(0)
        axE.text(mu[0] + off[0], mu[1] + off[1], lab, fontsize=5.8, color=col, ha="center", va="center",
                 zorder=10, path_effects=[pe.withStroke(linewidth=1.7, foreground="white")])
    axE.legend(loc="lower left", bbox_to_anchor=(-0.01, -0.01), ncol=2, frameon=False, markerscale=4.5,
               handletextpad=0.25, columnspacing=0.7, fontsize=5.8, borderpad=0.1, labelspacing=0.25)
    for a in (axD, axE):
        a.set_xlim(x0 - px, x1 + px); a.set_ylim(y0 - py, y1 + py)
        a.set_xticks([]); a.set_yticks([]); a.set_xlabel("UMAP 1"); a.set_ylabel("UMAP 2")
        a.set_box_aspect(1); a.spines["top"].set_visible(True); a.spines["right"].set_visible(True)
    axD.text(-0.16, 1.06, "D.", transform=axD.transAxes, fontsize=9, fontweight="bold", va="top")
    axE.text(-0.16, 1.06, "E.", transform=axE.transAxes, fontsize=9, fontweight="bold", va="top")
    bp = os.path.join(FIG, "Fig4DE_traj_pseudotime")
    for _e, _kw in (("_preview.png", dict(dpi=400)), (".pdf", dict()), (".svg", dict()),
                   (".tiff", dict(dpi=1200, pil_kwargs={"compression": "tiff_lzw"}))):
        fig.savefig(bp + _e, **{**_kw})
        print("  %-52s %6.2f MB" % (os.path.basename(bp + _e), os.path.getsize(bp + _e) / 1e6))
    plt.close(fig)
    print("saved", p)

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\traj_lines.log", "w", encoding="utf8", buffering=1)
    main()
