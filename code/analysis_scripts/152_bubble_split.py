# -*- coding: utf-8 -*-
"""152_bubble_split.py - 把注释证据气泡图拆成两张独立图：eLife 血细胞 / 肝胰腺 12 类"""
import io, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np, pandas as pd
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
OUT = r"D:\桌面\single cell\analysis_output"; FIG = os.path.join(OUT, "figures_pub")
N = pd.read_csv(os.path.join(OUT, "bubble2_n_supported.csv"), index_col=0)
A = pd.read_csv(os.path.join(OUT, "bubble2_auc.csv"), index_col=0)
ROWS = list(N.index)

def draw(cols, title, fname, figsize, rot=45, xfs=6.2, scale=3.6, legend=(5, 20, 50), lmax=None):
    n = N[cols].values; a = A[cols].values
    fig, ax = plt.subplots(figsize=figsize)
    for i in range(n.shape[0]):
        for j in range(n.shape[1]):
            if not np.isfinite(a[i, j]):
                continue
            ax.scatter(j, i, s=14 + n[i, j] * scale,
                       c=[plt.cm.viridis((np.clip(a[i, j], 0.5, 1.0) - 0.5) / 0.5)],
                       edgecolors="#555555", linewidths=0.3, zorder=3)
    ax.set_xlim(-0.6, len(cols) - 0.4); ax.set_ylim(len(ROWS) - 0.4, -0.6)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([c.replace("eLife ", "").replace("HP ", "") for c in cols], rotation=rot, ha="right", fontsize=xfs)
    ax.set_yticks(range(len(ROWS))); ax.set_yticklabels(ROWS, fontsize=7)
    ax.set_title(title, fontsize=8.5, pad=6)
    ax.grid(True, ls=(0, (1, 2)), lw=0.3, color="#E5E5E5", zorder=0)
    for s in ax.spines.values(): s.set_visible(False)
    # colour bar
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=Normalize(0.5, 1.0))
    cb = fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.22)
    cb.set_label("Set-level AUROC", fontsize=6.5); cb.ax.tick_params(labelsize=5.5); cb.outline.set_visible(False)
    # size legend
    ax.text(len(cols) + 0.12, -0.30, "genes", fontsize=5.8, va="center", clip_on=False)
    for k, v in enumerate(legend):
        ax.scatter(len(cols) + 0.10, 0.35 + k * 1.05, s=14 + v * scale, c="none", edgecolors="#666666",
                   linewidths=0.4, clip_on=False)
        ax.text(len(cols) + 0.42, 0.35 + k * 1.05, "%d" % v, fontsize=5.8, va="center", clip_on=False)
    p = os.path.join(FIG, fname)
    fig.savefig(p, dpi=400, bbox_inches="tight", pad_inches=0.03); plt.close(fig)
    print("saved", p)

def main():
    els = [c for c in N.columns if c.startswith("eLife")]
    hps = [c for c in N.columns if c.startswith("HP")]
    draw(els, "eLife hemocyte subtypes", "Fig_annotation_bubbles_blood_preview.png", (4.6, 3.4))
    draw(hps, "Hepatopancreas cell types (PMID 42004253; bubble size scaled to 1-4 markers)",
         "Fig_annotation_bubbles_hepatopancreas_preview.png", (6.6, 3.4), xfs=5.8, scale=30, legend=(1, 2, 4))

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\bubble_split.log", "w", encoding="utf8", buffering=1)
    main()
