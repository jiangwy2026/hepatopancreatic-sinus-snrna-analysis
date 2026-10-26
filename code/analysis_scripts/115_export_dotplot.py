# -*- coding: utf-8 -*-
"""115_export_dotplot.py - export the approved marker dot plot in all four formats."""
import sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd, scanpy as sc
import scipy.sparse as sp
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"; FIG = os.path.join(OUT, "figures_pub")
ROWS = ["PH", "GH", "MH1", "MH2"]

def main():
    A = sc.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))
    base = {"2": "PH", "1": "GH", "0": "MH1", "3": "MH2"}
    ct = A.obs["vendor_cluster"].astype(str).map(base)
    H = A[ct.isin(ROWS).values].copy()
    H.obs["ct"] = pd.Categorical(ct[ct.isin(ROWS)].values, categories=ROWS)
    H.X = H.layers["counts"].copy(); sc.pp.normalize_total(H, target_sum=1e4); sc.pp.log1p(H)
    gpos = {g: i for i, g in enumerate(H.var_names.astype(str))}
    P = pd.read_csv(os.path.join(OUT, "blood4_dotplot_genes.csv"))
    P["id"] = P["gid"].astype(str)
    X = H[:, [gpos[g] for g in P.id]].X
    X = np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X)
    pct = np.zeros((len(P), len(ROWS))); mean = np.zeros_like(pct)
    for j, g in enumerate(ROWS):
        m = (H.obs["ct"] == g).values
        pct[:, j] = 100 * (X[m] > 0).mean(0); mean[:, j] = X[m].mean(0)
    z = (mean - mean.mean(1, keepdims=True)) / (mean.std(1, keepdims=True) + 1e-9)
    blocks = [("Prohemocyte", P[P.grp == "PH"]), ("Granulocyte", P[P.grp == "GH"]),
              ("Monocytic hemocyte", P[P.grp.isin(["MH1", "MH2"])])]
    order = np.array([i for blk in blocks for i in blk[1].index])
    labels, pcto, zo = P.sym.values[order], pct[order], z[order]

    fig, ax = plt.subplots(figsize=(7.4, 2.5))
    for i in range(len(labels)):
        for j in range(len(ROWS)):
            if pcto[i, j] <= 0.5: continue
            ax.scatter(i, j, s=3 + (pcto[i, j] / 100) ** 1.25 * 62,
                       c=[plt.cm.RdBu_r((np.clip(zo[i, j], -2, 2) + 2) / 4)],
                       edgecolors="#555555", linewidths=0.25, zorder=3)
    ax.set_xlim(-0.7, len(labels) - 0.3); ax.set_ylim(len(ROWS) - 0.45, -0.55)
    ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=5.2, style="italic")
    ax.set_yticks(range(len(ROWS))); ax.set_yticklabels(ROWS, fontsize=7)
    for t in ax.get_yticklabels(): t.set_fontweight("bold")
    ax.tick_params(axis="both", length=0, pad=2)
    ax.grid(True, which="major", ls=(0, (1, 2)), lw=0.35, color="#CCCCCC", zorder=0)
    for s in ax.spines.values(): s.set_visible(False)
    k = 0
    for name, blk in blocks:
        a, b = k, k + len(blk) - 1; k = b + 1
        ax.plot([a - 0.45, b + 0.45], [-0.72, -0.72], color="black", lw=0.9, clip_on=False, zorder=4)
        ax.text((a + b) / 2, -1.05, name, ha="center", va="bottom", fontsize=6.4, weight="bold")
    pos = ax.get_position()
    cax = fig.add_axes([pos.x1 + 0.012, pos.y0 + 0.42 * pos.height, 0.010, 0.34 * pos.height])
    sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=plt.Normalize(-2, 2))
    cb = fig.colorbar(sm, cax=cax); cb.set_label("avg.exp.scale", fontsize=5.6, labelpad=1)
    cb.ax.tick_params(labelsize=5, length=1.2); cb.outline.set_visible(False)
    sx = pos.x1 + 0.055
    for kk, p in enumerate((20, 40, 60, 80)):
        y = 0.95 - kk * 0.16
        fig.text(sx, pos.y0 + y * pos.height, "●", fontsize=4 + (p / 100) ** 1.25 * 62 ** 0.5, va="center", color="#666666")
        fig.text(sx + 0.030, pos.y0 + y * pos.height, "%d" % p, fontsize=5.2, va="center")
    fig.text(sx + 0.010, pos.y0 + 1.02 * pos.height, "pct.exp", fontsize=5.6, va="bottom")

    bp = os.path.join(FIG, "Fig_blood4_dotplot_refstyle")
    for e, kw in (("_preview.png", dict(dpi=400)), (".pdf", {}), (".svg", {}),
                  (".tiff", dict(dpi=1200, pil_kwargs={"compression": "tiff_lzw"}))):
        fig.savefig(bp + e, bbox_inches="tight", pad_inches=0.03, **kw)
        print("  %-46s %6.2f MB" % (os.path.basename(bp + e), os.path.getsize(bp + e) / 1e6))
    plt.close(fig)

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\export_dotplot.log", "w", encoding="utf8", buffering=1)
    main()
