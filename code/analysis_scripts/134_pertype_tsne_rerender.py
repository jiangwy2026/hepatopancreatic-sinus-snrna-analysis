# -*- coding: utf-8 -*-
"""134_pertype_tsne_rerender.py - 用已保存的 t-SNE 重画（英文标签，避免 CJK 缺字形）"""
import io, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd, anndata as ad
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
OUT = r"D:\桌面\single cell\analysis_output"; FIG = os.path.join(OUT, "figures_pub")
ORD = ["PH", "GH", "MH1", "MH2", "Epi", "F", "R", "B", "NEU"]
COL = {"PH": "#B2182B", "GH": "#E08214", "MH1": "#FDB863", "MH2": "#8C510A", "Epi": "#7B3294",
       "F": "#006D2C", "R": "#41AB5D", "B": "#78C679", "NEU": "#2166AC"}

def main():
    B = ad.read_h5ad(os.path.join(OUT, "adata_pertypeQC_harmony_tsne.h5ad"))
    U = np.asarray(B.obsm["X_tsne"]); s = B.obs["sample"].astype(str).values
    c = B.obs["ct"].astype(str).values; added = B.obs["added"].values.astype(bool)
    thr = pd.read_csv(os.path.join(OUT, "pertype_thresholds_used.csv"))
    print(thr.to_string(index=False))

    # 主图：control vs treated
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.9))
    for ax, smp, title in zip(axes, ("blank", "treated"), ("Control (0 h)", "Hemospasia (6 h)")):
        m = s == smp
        for g in ORD:
            mm = m & (c == g)
            if mm.sum(): ax.scatter(U[mm, 0], U[mm, 1], s=1.2, c=COL[g], linewidths=0, rasterized=True,
                                    label="%s (%d)" % (g, int(mm.sum())))
        ax.set_title("%s   n = %s" % (title, format(int(m.sum()), ",")), fontsize=8)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("t-SNE 1"); ax.set_ylabel("t-SNE 2"); ax.set_box_aspect(1)
        for sp_ in ax.spines.values(): sp_.set_color("#999999")
    axes[1].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, markerscale=4,
                   fontsize=6.2, handletextpad=0.3, labelspacing=0.3)
    fig.suptitle("Per-cell-type QC (pct_mt < median + 3xMAD; n_genes >= 200) - %s cells" % format(len(s), ","), fontsize=8, y=1.0)
    p = os.path.join(FIG, "Fig_QC_pertype_tsne_blank_treated_preview.png")
    fig.savefig(p, dpi=400, bbox_inches="tight", pad_inches=0.03); plt.close(fig); print("saved", p)

    # 诊断：新纳入细胞落在哪里
    fig, axes = plt.subplots(1, 3, figsize=(11.6, 3.9))
    for ax, smp in zip(axes[:2], ("blank", "treated")):
        m = s == smp; old = m & ~added; new = m & added
        ax.scatter(U[old, 0], U[old, 1], s=1.2, c="#D9D9D9", linewidths=0, rasterized=True, label="kept by vendor QC (%d)" % old.sum())
        ax.scatter(U[new, 0], U[new, 1], s=1.4, c="#2166AC", linewidths=0, rasterized=True, label="added by per-type QC (%d)" % new.sum())
        ax.set_title("%s: %s newly added" % (smp, format(int(new.sum()), ",")), fontsize=8)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_box_aspect(1)
        ax.legend(fontsize=6, frameon=False, markerscale=4); ax.set_xlabel("t-SNE 1"); ax.set_ylabel("t-SNE 2")
        for sp_ in ax.spines.values(): sp_.set_color("#999999")
    ax = axes[2]
    ax.scatter(U[~added, 0], U[~added, 1], s=1.2, c="#D9D9D9", linewidths=0, rasterized=True, label="vendor-kept")
    mm = added & (B.obs.pct_counts_mt.values > 20)
    ax.scatter(U[mm, 0], U[mm, 1], s=1.4, c="#B2182B", linewidths=0, rasterized=True,
               label="newly added, pct_mt > 20%% (%s)" % format(int(mm.sum()), ","))
    ax.set_title("Newly added cells with high mitochondrial fraction", fontsize=8)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_box_aspect(1); ax.legend(fontsize=6, frameon=False, markerscale=4)
    ax.set_xlabel("t-SNE 1"); ax.set_ylabel("t-SNE 2")
    for sp_ in ax.spines.values(): sp_.set_color("#999999")
    p2 = os.path.join(FIG, "Fig_QC_pertype_added_cells_preview.png")
    fig.savefig(p2, dpi=400, bbox_inches="tight", pad_inches=0.03); plt.close(fig); print("saved", p2)

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\pertype_tsne.log", "w", encoding="utf8", buffering=1)
    main()
