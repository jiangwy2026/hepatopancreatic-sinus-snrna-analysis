# -*- coding: utf-8 -*-
"""
109_tsne_genes_and_E.py
(A) single-gene t-SNE feature plots over ALL cells (dark = high expression)
(B) the lineage figure (E) redrawn on the t-SNE embedding
"""
import sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Ellipse
from matplotlib.gridspec import GridSpec
import numpy as np, pandas as pd, scanpy as sc
from scipy.interpolate import UnivariateSpline
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"; FIG = os.path.join(OUT, "figures_pub")
BASE = {"2": "PH", "1": "GH", "0": "MH1", "3": "MH2", "4": "Epi", "5": "F", "6": "R", "7": "B",
        "8": "NE", "9": "Neu"}
ORD9 = ["PH", "GH", "MH1", "MH2", "Epi", "F", "R", "B", "NEU"]
BLOOD = ["PH", "GH", "MH1", "MH2"]
COL = {"PH": "#B2182B", "GH": "#E08214", "MH1": "#FDB863", "MH2": "#8C510A"}
GENES = [("GRIK2", "PH"), ("ccnd2", "PH"), ("Pxt", "GH"), ("FCN2", "GH"),
         ("PEN-3", "MH1"), ("THBS4", "MH1"), ("ITGA2B", "MH2"), ("SDC2", "MH2")]

def ell(pts, nsd):
    mu = pts.mean(0); cov = np.cov(pts.T)
    ev, evec = np.linalg.eigh(cov); o = ev.argsort()[::-1]; ev, evec = ev[o], evec[:, o]
    return mu, (2 * nsd * np.sqrt(ev[0]), 2 * nsd * np.sqrt(ev[1])), np.degrees(np.arctan2(evec[1, 0], evec[0, 0]))

def main():
    T = sc.read_h5ad(os.path.join(OUT, "adata_soupx_hvg_tsne.h5ad"))
    tsne = np.asarray(T.obsm["X_tsne"])
    A = sc.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))
    A = A[T.obs_names].copy()
    ct = A.obs["vendor_cluster"].astype(str).map(BASE).replace({"NE": "NEU", "Neu": "NEU"}).astype(str).values
    A.X = A.layers["counts"].copy(); sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
    sym2i = {}
    for i, s in enumerate(A.var["sym"].astype(str).values): sym2i.setdefault(s, i)
    miss = [g for g, _ in GENES if g not in sym2i]
    print("genes missing:", miss, flush=True)

    # ---------- A: single-gene feature plots on t-SNE
    n = len(GENES)
    fig, axes = plt.subplots(2, 4, figsize=(9.6, 4.9))
    for k, (g, grp) in enumerate(GENES):
        ax = axes.ravel()[k]
        if g not in sym2i:
            ax.axis("off"); continue
        v = np.asarray(A[:, sym2i[g]].X.todense()).ravel()
        o = np.argsort(-v)
        vmax = np.percentile(v[v > 0], 98) if (v > 0).any() else 1
        sc0 = ax.scatter(tsne[o, 0], tsne[o, 1], c=v[o], s=1.0, cmap="Reds", vmin=0, vmax=vmax,
                         linewidths=0, rasterized=True)
        for gg in BLOOD:
            m = ct == gg
            if m.sum() < 20: continue
            cx, cy = np.median(tsne[m, 0]), np.median(tsne[m, 1])
            ax.text(cx, cy, gg, fontsize=5.4, color="#404040", ha="center", va="center",
                    path_effects=[pe.withStroke(linewidth=1.3, foreground="white")], zorder=5)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_box_aspect(1)
        ax.set_title("%s  (%s)" % (g, grp), fontsize=7.2, pad=3)
        for s_ in ax.spines.values(): s_.set_color("#999999")
        cax = ax.inset_axes([0.86, 0.02, 0.11, 0.30])
        cb = fig.colorbar(sc0, cax=cax); cb.ax.tick_params(labelsize=4.2, length=1.0); cb.outline.set_visible(False)
    fig.suptitle("Single-gene expression on t-SNE (all %d cells; dark = high)" % T.n_obs, fontsize=8, y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    p1 = os.path.join(FIG, "Fig_tsne_genefeatures_preview.png")
    fig.savefig(p1, dpi=400, bbox_inches="tight", pad_inches=0.03); plt.close(fig)
    print("saved", p1, flush=True)

    # ---------- B: lineage figure (E) on t-SNE
    m4 = np.isin(ct, BLOOD)
    U = tsne[m4]; c4 = ct[m4]
    cent = {g: U[c4 == g].mean(0) for g in BLOOD}
    x0, x1 = np.percentile(U[:, 0], [0.3, 99.7]); y0, y1 = np.percentile(U[:, 1], [0.3, 99.7])
    px = 0.04 * (x1 - x0); py = 0.04 * (y1 - y0)
    fig, ax = plt.subplots(figsize=(3.5, 3.4))
    for g in BLOOD:
        m = c4 == g
        ax.scatter(U[m, 0], U[m, 1], s=1.3, c=COL[g], linewidths=0, rasterized=True,
                   label="%s (%d)" % (g, int(m.sum())))
    for m, off in [((c4 == "GH"), (0.4, 4.3)), ((c4 == "MH1") | (c4 == "MH2"), (-5.4, -5.9)), ((c4 == "PH"), (1.4, 4.6))]:
        mu, (w, h), ang = ell(U[m], 1.7)
        ax.add_patch(Ellipse(mu, w, h, angle=ang, fill=False, ls=(0, (3.0, 2.4)), lw=0.55, ec="#9A9A9A", zorder=5))
    for path, c in [(("PH", "GH"), "#2E6DA4"), (("PH", "MH2", "MH1"), "#7F4A12")]:
        p = np.vstack([cent[g] for g in path])
        t = np.arange(len(p), dtype=float)
        fx = UnivariateSpline(t, p[:, 0], k=min(3, len(p) - 1)); fy = UnivariateSpline(t, p[:, 1], k=min(3, len(p) - 1))
        tt = np.linspace(0, len(p) - 1, 220); curve = np.c_[fx(tt), fy(tt)]
        ax.plot(curve[:, 0], curve[:, 1], lw=1.5, color=c, solid_capstyle="round", zorder=8,
                path_effects=[pe.withStroke(linewidth=2.7, foreground="white")])
        for i in (len(curve) // 2, len(curve) - 6):
            if i + 4 >= len(curve): i = len(curve) - 8
            ax.annotate("", xy=curve[i + 4], xytext=curve[i],
                        arrowprops=dict(arrowstyle="-|>,head_length=0.55,head_width=0.24", color=c, lw=0,
                                        shrinkA=0, shrinkB=0, mutation_scale=11), zorder=9)
    for key, lab, col, off in [("GH", "Granulocyte lineage", "#2E6DA4", (0.4, 4.3)),
                               ("MH", "Monocytic hemocyte lineage", "#7F4A12", (-5.4, -5.9)),
                               ("PH", "Prohemocyte", "#B2182B", (1.4, 4.6))]:
        m = (c4 == "GH") if key == "GH" else ((c4 == "MH1") | (c4 == "MH2")) if key == "MH" else (c4 == "PH")
        mu = U[m].mean(0)
        ax.text(mu[0] + off[0], mu[1] + off[1], lab, fontsize=5.9, color=col, ha="center", va="center",
                zorder=10, path_effects=[pe.withStroke(linewidth=1.7, foreground="white")])
    ax.legend(loc="lower left", bbox_to_anchor=(-0.01, -0.01), ncol=2, frameon=False, markerscale=4.5,
              handletextpad=0.25, columnspacing=0.7, fontsize=6, borderpad=0.1, labelspacing=0.25)
    ax.set_xlim(x0 - px, x1 + px); ax.set_ylim(y0 - py, y1 + py)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("t-SNE 1"); ax.set_ylabel("t-SNE 2")
    ax.set_box_aspect(1); ax.spines["top"].set_visible(True); ax.spines["right"].set_visible(True)
    p2 = os.path.join(FIG, "Fig4E_lineages_tsne_preview.png")
    fig.savefig(p2, dpi=400, bbox_inches="tight", pad_inches=0.03); plt.close(fig)
    print("saved", p2, flush=True)

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\tsne_figs.log", "w", encoding="utf8", buffering=1)
    main()
