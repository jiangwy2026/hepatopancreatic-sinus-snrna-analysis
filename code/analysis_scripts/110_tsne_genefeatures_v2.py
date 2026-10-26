# -*- coding: utf-8 -*-
"""110_tsne_genefeatures_v2.py - single-gene t-SNE feature plots with grey background and
coloured expressing cells; each panel also reports the %% of expressing cells."""
import sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np, pandas as pd, scanpy as sc
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"; FIG = os.path.join(OUT, "figures_pub")
BASE = {"2": "PH", "1": "GH", "0": "MH1", "3": "MH2", "4": "Epi", "5": "F", "6": "R", "7": "B", "8": "NE", "9": "Neu"}
BLOOD = ["PH", "GH", "MH1", "MH2"]
GENES = [("GRIK2", "PH"), ("ccnd2", "PH"), ("Pxt", "GH"), ("A2M", "GH"),
         ("PEN-3", "MH1"), ("THBS4", "MH1"), ("ITGA2B", "MH2"), ("SDC2", "MH2")]

def main():
    T = sc.read_h5ad(os.path.join(OUT, "adata_soupx_hvg_tsne.h5ad"))
    tsne = np.asarray(T.obsm["X_tsne"])
    A = sc.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))[T.obs_names].copy()
    ct = A.obs["vendor_cluster"].astype(str).map(BASE).replace({"NE": "NEU", "Neu": "NEU"}).astype(str).values
    A.X = A.layers["counts"].copy(); sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
    sym2i = {}
    for i, s in enumerate(A.var["sym"].astype(str).values): sym2i.setdefault(s, i)

    rows = []
    for g, grp in GENES:
        v = np.asarray(A[:, sym2i[g]].X.todense()).ravel()
        rows.append({"gene": g, "group": grp, "pct_all": round(100 * (v > 0).mean(), 1),
                     **{gg: round(100 * (v[ct == gg] > 0).mean(), 1) for gg in BLOOD}})
    tab = pd.DataFrame(rows)
    print(tab.to_string(index=False)); tab.to_csv(os.path.join(OUT, "tsne_feature_gene_pct.csv"), index=False)

    fig, axes = plt.subplots(2, 4, figsize=(9.6, 5.1))
    for k, (g, grp) in enumerate(GENES):
        ax = axes.ravel()[k]
        v = np.asarray(A[:, sym2i[g]].X.todense()).ravel()
        ax.scatter(tsne[:, 0], tsne[:, 1], c="#E8E8E8", s=0.8, linewidths=0, rasterized=True)
        pos = v > 0
        vmax = np.percentile(v[pos], 95) if pos.any() else 1.0
        o = np.argsort(-v)
        sel = np.array([i for i in o if v[i] > 0])
        sc0 = ax.scatter(tsne[sel, 0], tsne[sel, 1], c=v[sel], cmap="Reds", vmin=0, vmax=vmax,
                         s=1.8, linewidths=0, rasterized=True)
        for gg in BLOOD:
            m = ct == gg
            cx, cy = np.median(tsne[m, 0]), np.median(tsne[m, 1])
            ax.text(cx, cy, gg, fontsize=5.4, color="#3A3A3A", ha="center", va="center", zorder=5,
                    path_effects=[pe.withStroke(linewidth=1.4, foreground="white")])
        ax.set_xticks([]); ax.set_yticks([]); ax.set_box_aspect(1)
        ax.set_title("%s  (%s)   %0.1f%% cells" % (g, grp, 100 * pos.mean()), fontsize=7.0, pad=3)
        for s_ in ax.spines.values(): s_.set_color("#999999")
        cax = ax.inset_axes([0.87, 0.03, 0.10, 0.26])
        cb = fig.colorbar(sc0, cax=cax); cb.ax.tick_params(labelsize=4.2, length=1.0); cb.outline.set_visible(False)
    fig.suptitle("Expression of single marker genes across all %d cells (t-SNE; grey = not detected, dark red = high)" % T.n_obs,
                 fontsize=8, y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    p = os.path.join(FIG, "Fig_tsne_genefeatures_preview.png")
    fig.savefig(p, dpi=400, bbox_inches="tight", pad_inches=0.03); plt.close(fig)
    print("saved", p)

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\tsne_figs.log", "w", encoding="utf8", buffering=1)
    main()
