# -*- coding: utf-8 -*-
"""190_violin_all_volcano_genes.py
All genes labelled in the early-PH volcano, shown together as violins:
positive: hemocyte TGM, ccnd2, SSPO, Notch-like;
negative: THBS4, PEN-3, crustin, WAP.
"""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
from scipy.stats import mannwhitneyu
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"
PREV = os.path.join(OUT, "figures_pub", "sc_summary_new_names", "individual_previews")
os.makedirs(PREV, exist_ok=True)
CTRL = "#4C78A8"; TRT = "#D6604D"
GENES = [
    ("ncbi_113823945", "hemocyte TGM", "pos"),
    ("ncbi_113814652", "ccnd2", "pos"),
    ("ncbi_113824171", "SSPO", "pos"),
    ("ncbi_113808957", "Notch-like", "pos"),
    ("ncbi_113824160", "THBS4", "neg"),
    ("ncbi_113808997", "PEN-3", "neg"),
    ("ncbi_113801825", "crustin", "neg"),
    ("ncbi_113810905", "WAP", "neg"),
]

def main():
    A = ad.read_h5ad(os.path.join(OUT, "adata_C2_within.h5ad"))
    A.X = A.layers["counts"].copy(); sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
    A.obs["pt_bin3"] = pd.qcut(A.obs["dpt_pseudotime"].astype(float), 3, labels=["P1", "P2", "P3"])
    E = A[A.obs["pt_bin3"] == "P1"].copy()
    gid2i = {str(g): i for i, g in enumerate(A.var_names.astype(str))}
    grp = E.obs["group"].astype(str).values
    rng = np.random.default_rng(3)

    fig = plt.figure(figsize=(7.2, 4.9), facecolor="white")
    gs = fig.add_gridspec(2, 4, left=0.10, right=0.985, top=0.86, bottom=0.15,
                          wspace=0.62, hspace=0.66)
    axes = []
    for k, (gid, label, kind) in enumerate(GENES):
        ax = fig.add_subplot(gs[k // 4, k % 4]); axes.append(ax)
        v = E[:, gid2i[gid]].X.toarray().ravel().astype(float)
        data = [v[grp == "blank"], v[grp == "treated"]]
        vp = ax.violinplot(data, positions=[1, 2], widths=0.62, showextrema=False, showmeans=False)
        for body, col in zip(vp["bodies"], (CTRL, TRT)):
            body.set_facecolor(col); body.set_edgecolor("none"); body.set_alpha(0.55)
        for pos, val, col in zip((1, 2), data, (CTRL, TRT)):
            j = rng.uniform(-0.07, 0.07, size=len(val))
            ax.scatter(pos + j, val, s=0.9, c=col, alpha=0.42, linewidths=0, rasterized=True)
            ax.scatter([pos], [np.median(val)], s=11, c="white", edgecolors=col, linewidths=0.6, zorder=5)
        p = mannwhitneyu(data[1], data[0], alternative="two-sided").pvalue
        ymax = max(np.max(data[0]), np.max(data[1]))
        top = ymax + (0.25 if ymax > 0 else 0.1)
        ax.set_ylim(-0.08, top * 1.08 + 0.05)
        ax.plot([1, 2], [top, top], color="#333333", lw=0.45)
        ax.text(1.5, top + 0.025 * top, "P = %.1e" % p, ha="center", va="bottom", fontsize=4.9)
        ax.set_xticks([1, 2]); ax.set_xticklabels(["Control", "Hemospasia"], fontsize=4.7, rotation=20, ha="right")
        ax.tick_params(axis="y", labelsize=4.8, length=1.1, width=0.45)
        ax.tick_params(axis="x", length=0, pad=1)
        color = "#B2182B" if kind == "pos" else "#2166AC"
        ax.set_title(label, fontsize=5.8, color=color, pad=4)
        if k % 4 == 0:
            ax.set_ylabel("Log expression", fontsize=5.4)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.text(-0.28, 1.04, "abcdefgh"[k], transform=ax.transAxes, fontsize=7.2,
                fontweight="bold", va="top")

    fig.text(0.025, 0.69, "Activated in treated PH", rotation=90, ha="center", va="center",
             fontsize=5.8, color="#B2182B", weight="bold")
    fig.text(0.025, 0.30, "Suppressed in treated PH", rotation=90, ha="center", va="center",
             fontsize=5.8, color="#2166AC", weight="bold")
    fig.legend(handles=[Patch(color=CTRL, label="Control"), Patch(color=TRT, label="Hemospasia")],
               loc="upper center", bbox_to_anchor=(0.54, 0.925), ncol=2, frameon=False,
               fontsize=5.4, handlelength=0.9, handletextpad=0.3, columnspacing=1.1)
    fig.suptitle("Early PH: genes labelled in the volcano", fontsize=7.0, x=0.54, y=0.965)
    p = os.path.join(PREV, "04c_violin_all_volcano_genes_preview.png")
    fig.savefig(p, dpi=400, bbox_inches="tight", pad_inches=0.04, facecolor="white")
    print("saved", p)
    plt.close(fig)

if __name__ == "__main__":
    main()
