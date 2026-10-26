# -*- coding: utf-8 -*-
"""102_figE_only.py - single-panel (E) figure: functional groups + lineage ellipses + two trajectory lines."""
import os, sys, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Ellipse
import numpy as np
from scipy.interpolate import UnivariateSpline
import anndata as ad

OUT = r"D:\桌面\single cell\analysis_output"; FIG = os.path.join(OUT, "figures_pub")
sys.path.insert(0, OUT)
from pub_style import apply_pub_style
apply_pub_style(7)
ORD = ["PH", "GH", "MH1", "MH2"]
COL = {"PH": "#B2182B", "GH": "#E08214", "MH1": "#FDB863", "MH2": "#8C510A"}
LIN = [("PH", "GH", "#2E6DA4"), ("PH", "MH2", "MH1", "#7F4A12")]     # granulocyte / monocytic

def ell(pts, nsd):
    mu = pts.mean(0); cov = np.cov(pts.T)
    ev, evec = np.linalg.eigh(cov); o = ev.argsort()[::-1]; ev, evec = ev[o], evec[:, o]
    return mu, (2 * nsd * np.sqrt(ev[0]), 2 * nsd * np.sqrt(ev[1])), np.degrees(np.arctan2(evec[1, 0], evec[0, 0]))

def main():
    A = ad.read_h5ad(os.path.join(OUT, "adata_hemocyte_traj_geodesic.h5ad"))
    ct = A.obs["ct"].astype(str).values
    U = np.asarray(A.obsm["X_umap"])
    cent = {g: U[ct == g].mean(0) for g in ORD}
    x0, x1 = np.percentile(U[:, 0], [0.3, 99.7]); y0, y1 = np.percentile(U[:, 1], [0.3, 99.7])
    px = 0.04 * (x1 - x0); py = 0.04 * (y1 - y0)

    fig, ax = plt.subplots(figsize=(3.5, 3.4))
    for g in ORD:
        m = ct == g
        ax.scatter(U[m, 0], U[m, 1], s=1.2, c=COL[g], linewidths=0, rasterized=True,
                   label="%s (%d)" % (g, int(m.sum())))
    for m, off in [((ct == "GH"), (0.4, 4.3)), ((ct == "MH1") | (ct == "MH2"), (-5.4, -5.9)), ((ct == "PH"), (1.4, 4.6))]:
        mu, (w, h), ang = ell(U[m], 1.7)
        ax.add_patch(Ellipse(mu, w, h, angle=ang, fill=False, ls=(0, (3.0, 2.4)), lw=0.55, ec="#9A9A9A", zorder=5))
    for spec in LIN:
        path = spec[:-1]; c = spec[-1]
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
        m = (ct == "GH") if key == "GH" else ((ct == "MH1") | (ct == "MH2")) if key == "MH" else (ct == "PH")
        mu = U[m].mean(0)
        ax.text(mu[0] + off[0], mu[1] + off[1], lab, fontsize=5.9, color=col, ha="center", va="center",
                zorder=10, path_effects=[pe.withStroke(linewidth=1.7, foreground="white")])
    ax.legend(loc="lower left", bbox_to_anchor=(-0.01, -0.01), ncol=2, frameon=False, markerscale=4.5,
              handletextpad=0.25, columnspacing=0.7, fontsize=6, borderpad=0.1, labelspacing=0.25)
    ax.set_xlim(x0 - px, x1 + px); ax.set_ylim(y0 - py, y1 + py)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    ax.set_box_aspect(1); ax.spines["top"].set_visible(True); ax.spines["right"].set_visible(True)
    bp = os.path.join(FIG, "Fig4E_lineages")
    for _e, _kw in (("_preview.png", dict(dpi=400)), (".pdf", {}), (".svg", {}),
                    (".tiff", dict(dpi=1200, pil_kwargs={"compression": "tiff_lzw"}))):
        fig.savefig(bp + _e, bbox_inches="tight", **_kw)
        print("  %-52s %6.2f MB" % (os.path.basename(bp + _e), os.path.getsize(bp + _e) / 1e6))
    plt.close(fig)
    print("saved", p)

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\figE.log", "w", encoding="utf8", buffering=1)
    main()
