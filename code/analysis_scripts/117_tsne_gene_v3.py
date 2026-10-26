# -*- coding: utf-8 -*-
"""117_tsne_gene_v3.py - single-gene t-SNE panels: grey background, green expression,
ONE shared colour bar placed outside the panels so no cells are hidden."""
import sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import ListedColormap
import numpy as np, pandas as pd, scanpy as sc
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"; FIG = os.path.join(OUT, "figures_pub")
BASE = {"2": "PH", "1": "GH", "0": "MH1", "3": "MH2", "4": "Epi", "5": "F", "6": "R", "7": "B", "8": "NE", "9": "Neu"}
ROWS = ["PH", "GH", "MH1", "MH2"]
GENES = [("PPAF2", "ncbi_113812814", "PH"), ("Grik2", "ncbi_113823802", "PH"),
         ("A2M", "ncbi_113809687", "GH"), ("Pxt", "ncbi_113820150", "GH"),
         ("LOC113819543", "ncbi_113819543", "MH1"), ("THBS4", "ncbi_113824160", "MH1"),
         ("LOC113810024", "ncbi_113810024", "MH2"), ("ITGA2B", "ncbi_113815644", "MH2")]
GREY = "#C9C9C9"
GCMAP = ListedColormap(plt.cm.Greens(np.linspace(0.28, 0.95, 256)))

def main():
    T = sc.read_h5ad(os.path.join(OUT, "adata_soupx_hvg_tsne.h5ad"))
    tsne = np.asarray(T.obsm["X_tsne"])
    B = sc.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))[T.obs_names].copy()
    ct = B.obs["vendor_cluster"].astype(str).map(BASE).replace({"NE": "NEU", "Neu": "NEU"}).astype(str).values
    B.X = B.layers["counts"].copy(); sc.pp.normalize_total(B, target_sum=1e4); sc.pp.log1p(B)
    bid = {g: i for i, g in enumerate(B.var_names.astype(str))}

    fig, axes = plt.subplots(2, 4, figsize=(9.6, 5.0))
    im_for_cb = None
    for k, (g, gid, grp) in enumerate(GENES):
        ax = axes.ravel()[k]
        v = np.asarray(B[:, bid[gid]].X.todense()).ravel()
        ax.scatter(tsne[:, 0], tsne[:, 1], c=GREY, s=0.9, linewidths=0, rasterized=True)
        pos = np.where(v > 0)[0]
        if len(pos):
            vmax = np.percentile(v[pos], 95)
            rel = np.clip(v[pos] / vmax, 0, 1)
            o = np.argsort(rel)
            im_for_cb = ax.scatter(tsne[pos[o], 0], tsne[pos[o], 1], c=rel[o], cmap=GCMAP,
                                   vmin=0, vmax=1, s=2.0, linewidths=0, rasterized=True)
        for gg in ROWS:
            m = ct == gg
            ax.text(np.median(tsne[m, 0]), np.median(tsne[m, 1]), gg, fontsize=5.4, color="#2B2B2B",
                    ha="center", va="center", zorder=6,
                    path_effects=[pe.withStroke(linewidth=1.6, foreground="white")])
        ax.set_xticks([]); ax.set_yticks([]); ax.set_box_aspect(1)
        ax.set_title("%s  (%s)   %0.1f%% cells" % (g, grp, 100 * (v > 0).mean()), fontsize=7.0, pad=3)
        for s_ in ax.spines.values(): s_.set_color("#999999")
    fig.tight_layout(rect=(0, 0.05, 1, 0.975))
    cax = fig.add_axes([0.72, 0.035, 0.18, 0.016])
    cb = fig.colorbar(im_for_cb, cax=cax, orientation="horizontal")
    cb.set_label("Relative expression (per gene, 0-1)", fontsize=6.5, labelpad=2)
    cb.ax.tick_params(labelsize=5.5, length=1.2); cb.outline.set_visible(False)
    fig.suptitle("Single-gene expression across all %d cells (t-SNE; grey = not detected, dark green = high)" % T.n_obs,
                 fontsize=8, y=0.995)
    p = os.path.join(FIG, "Fig_tsne_genefeatures_preview.png")
    fig.savefig(p, dpi=400, bbox_inches="tight", pad_inches=0.03); plt.close(fig)
    print("saved", p)

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\tsne_gene_v3.log", "w", encoding="utf8", buffering=1)
    main()
