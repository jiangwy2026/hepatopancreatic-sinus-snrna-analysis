# -*- coding: utf-8 -*-
"""118_figE_tsne_allgroups.py - lineage figure on the t-SNE keeping ALL subpopulations:
the five non-blood groups in grey, the four blood groups coloured, with lineage ellipses and curves."""
import sys, os, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Ellipse
import numpy as np, pandas as pd, scanpy as sc
from scipy.interpolate import UnivariateSpline
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"; FIG = os.path.join(OUT, "figures_pub")
BASE = {"2": "PH", "1": "GH", "0": "MH1", "3": "MH2", "4": "Epi", "5": "F", "6": "R", "7": "B", "8": "NE", "9": "Neu"}
BLOOD = ["PH", "GH", "MH1", "MH2"]
COL = {"PH": "#B2182B", "GH": "#E08214", "MH1": "#FDB863", "MH2": "#8C510A"}
GREY = "#CFCFCF"

def ell(pts, nsd):
    mu = pts.mean(0); cov = np.cov(pts.T)
    ev, evec = np.linalg.eigh(cov); o = ev.argsort()[::-1]; ev, evec = ev[o], evec[:, o]
    return mu, (2 * nsd * np.sqrt(ev[0]), 2 * nsd * np.sqrt(ev[1])), np.degrees(np.arctan2(evec[1, 0], evec[0, 0]))

def main():
    T = sc.read_h5ad(os.path.join(OUT, "adata_soupx_hvg_tsne.h5ad"))
    U = np.asarray(T.obsm["X_tsne"])
    ct = T.obs["vendor_cluster"].astype(str).map(BASE).replace({"NE": "NEU", "Neu": "NEU"}).astype(str).values
    other = ~np.isin(ct, BLOOD)
    print("all cells:", len(ct), "| blood:", int((~other).sum()), "| other:", int(other.sum()),
          "|", pd.Series(ct[other]).value_counts().to_dict(), flush=True)

    fig, ax = plt.subplots(figsize=(3.6, 3.5))
    ax.scatter(U[other, 0], U[other, 1], s=1.0, c=GREY, linewidths=0, rasterized=True,
               label="Other subpopulations (%d)" % int(other.sum()), zorder=2)
    for g in BLOOD:
        m = ct == g
        ax.scatter(U[m, 0], U[m, 1], s=1.3, c=COL[g], linewidths=0, rasterized=True,
                   label="%s (%d)" % (g, int(m.sum())), zorder=3)
    for m, in [((ct == "GH"),), ((ct == "MH1") | (ct == "MH2"),), ((ct == "PH"),)]:
        mu, (w, h), ang = ell(U[m], 1.7)
        ax.add_patch(Ellipse(mu, w, h, angle=ang, fill=False, ls=(0, (3.0, 2.4)), lw=0.55, ec="#9A9A9A", zorder=5))
    cent = {g: U[ct == g].mean(0) for g in BLOOD}
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
        m = (ct == "GH") if key == "GH" else ((ct == "MH1") | (ct == "MH2")) if key == "MH" else (ct == "PH")
        mu = U[m].mean(0)
        ax.text(mu[0] + off[0], mu[1] + off[1], lab, fontsize=5.9, color=col, ha="center", va="center",
                zorder=10, path_effects=[pe.withStroke(linewidth=1.7, foreground="white")])
    ax.legend(loc="lower left", bbox_to_anchor=(-0.01, -0.01), ncol=2, frameon=False, markerscale=4.5,
              handletextpad=0.25, columnspacing=0.7, fontsize=5.8, borderpad=0.1, labelspacing=0.22)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("t-SNE 1"); ax.set_ylabel("t-SNE 2")
    ax.set_box_aspect(1); ax.spines["top"].set_visible(True); ax.spines["right"].set_visible(True)
    p = os.path.join(FIG, "Fig4E_lineages_tsne_preview.png")
    fig.savefig(p, dpi=400, bbox_inches="tight", pad_inches=0.03); plt.close(fig)
    print("saved", p)

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\figE_tsne_all.log", "w", encoding="utf8", buffering=1)
    main()
