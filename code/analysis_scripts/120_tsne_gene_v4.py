# -*- coding: utf-8 -*-
"""120_tsne_gene_v4.py - single-gene t-SNE panels: per-panel colour bar on the right (outside the
axes so no cells are hidden) and more vertical space so titles never touch the panel above."""
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

    fig, axes = plt.subplots(2, 4, figsize=(10.4, 5.6),
                             gridspec_kw=dict(wspace=0.34, hspace=0.36, left=0.03, right=0.97, top=0.90, bottom=0.04))
    for k, (g, gid, grp) in enumerate(GENES):
        ax = axes.ravel()[k]
        v = np.asarray(B[:, bid[gid]].X.todense()).ravel()
        ax.scatter(tsne[:, 0], tsne[:, 1], c=GREY, s=0.9, linewidths=0, rasterized=True)
        pos = np.where(v > 0)[0]
        if len(pos):
            vmax = np.percentile(v[pos], 95)
            rel = np.clip(v[pos] / vmax, 0, 1)
            o = np.argsort(rel)
            sc0 = ax.scatter(tsne[pos[o], 0], tsne[pos[o], 1], c=rel[o], cmap=GCMAP, vmin=0, vmax=1,
                             s=2.0, linewidths=0, rasterized=True)
            cax = ax.inset_axes([1.035, 0.28, 0.042, 0.44])
            cb = fig.colorbar(sc0, cax=cax, orientation="vertical")
            cb.set_ticks([0, 0.5, 1.0])
            cb.ax.tick_params(labelsize=4.6, length=1.1, pad=1)
            cb.outline.set_visible(False)
            if k == 0:
                cb.set_label("Relative expression", fontsize=5.4, labelpad=2)
        for gg in ROWS:
            m = ct == gg
            ax.text(np.median(tsne[m, 0]), np.median(tsne[m, 1]), gg, fontsize=5.4, color="#2B2B2B",
                    ha="center", va="center", zorder=6,
                    path_effects=[pe.withStroke(linewidth=1.6, foreground="white")])
        ax.set_xticks([]); ax.set_yticks([]); ax.set_box_aspect(1)
        ax.set_title("%s  (%s)   %0.1f%% cells" % (g, grp, 100 * (v > 0).mean()), fontsize=7.2, pad=5)
        for s_ in ax.spines.values(): s_.set_color("#B0B0B0")
    fig.suptitle("Single-gene expression across all %d cells (t-SNE; grey = not detected, dark green = high)"
                 % T.n_obs, fontsize=8.5, y=0.985)
    bp = os.path.join(FIG, "Fig_tsne_genefeatures")
    for _e, _kw in (("_preview.png", dict(dpi=400)), (".pdf", dict(bbox_inches="tight", pad_inches=0.04)), (".svg", dict(bbox_inches="tight", pad_inches=0.04)),
                   (".tiff", dict(bbox_inches="tight", pad_inches=0.04, dpi=1200, pil_kwargs={"compression": "tiff_lzw"}))):
        fig.savefig(bp + _e, **_kw)
        print("  %-52s %6.2f MB" % (os.path.basename(bp + _e), os.path.getsize(bp + _e) / 1e6))
    plt.close(fig)
    print("saved", bp + "_preview.png")

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\tsne_gene_v4.log", "w", encoding="utf8", buffering=1)
    main()
