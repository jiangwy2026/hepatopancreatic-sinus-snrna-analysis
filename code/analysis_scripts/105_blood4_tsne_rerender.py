# -*- coding: utf-8 -*-
"""105_blood4_tsne_rerender.py - re-render the blood-lineage t-SNE with the legend below the panel."""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np
import anndata as ad
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
OUT = r"D:\桌面\single cell\analysis_output"; FIG = os.path.join(OUT, "figures_pub")
ORD = ["PH", "GH", "MH1", "MH2"]
COL = {"PH": "#B2182B", "GH": "#E08214", "MH1": "#FDB863", "MH2": "#8C510A"}

def main():
    A = ad.read_h5ad(os.path.join(OUT, "adata_blood4_tsne.h5ad"))
    ct = A.obs["ct"].astype(str).values
    U = np.asarray(A.obsm["X_tsne_blood"])
    fig, ax = plt.subplots(figsize=(3.6, 3.9))
    for g in ORD:
        m = ct == g
        ax.scatter(U[m, 0], U[m, 1], s=1.3, c=COL[g], linewidths=0, rasterized=True,
                   label="%s (%d)" % (g, int(m.sum())))
    for g in ORD:
        m = ct == g
        cx, cy = np.median(U[m, 0]), np.median(U[m, 1])
        ax.text(cx, cy, g, fontsize=7.5, weight="bold", ha="center", va="center", color="black",
                path_effects=[pe.withStroke(linewidth=2.0, foreground="white")])
    ax.set_xlabel("t-SNE 1"); ax.set_ylabel("t-SNE 2")
    ax.set_xticks([]); ax.set_yticks([]); ax.set_box_aspect(1)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.07), ncol=4, frameon=False, markerscale=5,
              handletextpad=0.25, columnspacing=0.8, fontsize=6, labelspacing=0.2)
    fig.tight_layout()
    p = os.path.join(FIG, "Fig_blood4_tsne_preview.png")
    fig.savefig(p, dpi=400, bbox_inches="tight", pad_inches=0.03); plt.close(fig)
    print("saved", p)

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\blood4_fig.log", "a", encoding="utf8", buffering=1)
    main()
