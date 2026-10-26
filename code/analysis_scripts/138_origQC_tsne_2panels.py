# -*- coding: utf-8 -*-
"""
138_origQC_tsne_2panels.py
用原阈值（厂商 QC, 9,658 细胞）出 Control / Treated 两张 t-SNE，
并做样本叠加图 + 组成柱状图 + 局部混合度，解释两组差异来源。
"""
import io, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np, pandas as pd, anndata as ad, scanpy as sc
from sklearn.neighbors import NearestNeighbors
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"; FIG = os.path.join(OUT, "figures_pub")
ORD = ["PH", "GH", "MH1", "MH2", "Epi", "F", "R", "B", "NEU"]
COL = {"PH": "#B2182B", "GH": "#E08214", "MH1": "#FDB863", "MH2": "#8C510A", "Epi": "#7B3294",
       "F": "#006D2C", "R": "#41AB5D", "B": "#78C679", "NEU": "#2166AC"}
BASE = {"2": "PH", "1": "GH", "0": "MH1", "3": "MH2", "4": "Epi", "5": "F", "6": "R", "7": "B", "8": "NEU", "9": "NEU"}

def main():
    A = ad.read_h5ad(os.path.join(OUT, "adata_soupx_hvg_tsne.h5ad"))
    ct = A.obs["vendor_cluster"].astype(str).map(BASE).values
    U = np.asarray(A.obsm["X_tsne"]); s = A.obs["sample"].astype(str).values
    print("cells:", len(s), "| blank", (s == "blank").sum(), "| treated", (s == "treated").sum())

    # A/B: 两张按类型着色
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.9))
    for ax, smp, title in zip(axes, ("blank", "treated"), ("Control (0 h)", "Hemospasia (6 h)")):
        m = s == smp
        for g in ORD:
            mm = m & (ct == g)
            if mm.sum():
                ax.scatter(U[mm, 0], U[mm, 1], s=1.2, c=COL[g], linewidths=0, rasterized=True,
                           label="%s (%d)" % (g, int(mm.sum())))
        ax.set_title("%s   n = %s" % (title, format(int(m.sum()), ",")), fontsize=8)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("t-SNE 1"); ax.set_ylabel("t-SNE 2"); ax.set_box_aspect(1)
        for sp_ in ax.spines.values(): sp_.set_color("#999999")
    axes[1].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, markerscale=4,
                   fontsize=6.2, handletextpad=0.3, labelspacing=0.3)
    fig.suptitle("Original QC (vendor: pct_mt < 10%%) - %s cells" % format(len(s), ","), fontsize=8, y=1.0)
    p = os.path.join(FIG, "Fig_QC_orig_tsne_blank_treated_preview.png")
    fig.savefig(p, dpi=400, bbox_inches="tight", pad_inches=0.03); plt.close(fig); print("saved", p)

    # C/D: 样本叠加 + 组成
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.9))
    ax = axes[0]
    for smp, c, lab in (("blank", "#2166AC", "Control (0 h)"), ("treated", "#B2182B", "Hemospasia (6 h)")):
        m = s == smp
        ax.scatter(U[m, 0], U[m, 1], s=0.9, c=c, linewidths=0, rasterized=True, alpha=0.45, label="%s (%d)" % (lab, m.sum()))
    ax.set_title("Both samples overlaid", fontsize=8)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("t-SNE 1"); ax.set_ylabel("t-SNE 2"); ax.set_box_aspect(1)
    ax.legend(fontsize=6, frameon=False, markerscale=5, loc="lower left")
    for sp_ in ax.spines.values(): sp_.set_color("#999999")
    ax = axes[1]
    t = pd.crosstab(ct, s).reindex(ORD)
    pct = t / t.sum(axis=0) * 100
    x = np.arange(len(ORD))
    ax.bar(x - 0.2, pct["blank"], 0.4, color="#2166AC", label="Control (0 h)")
    ax.bar(x + 0.2, pct["treated"], 0.4, color="#B2182B", label="Hemospasia (6 h)")
    for i, g in enumerate(ORD):
        ax.text(i, max(pct.loc[g, "blank"], pct.loc[g, "treated"]) + 1.2,
                "%.1fx" % (pct.loc[g, "treated"] / max(pct.loc[g, "blank"], 1e-9)), ha="center", fontsize=5.5)
    ax.set_xticks(x); ax.set_xticklabels(ORD, rotation=45, ha="right", fontsize=6.5)
    ax.set_ylabel("Cells (%) of sample", fontsize=7); ax.set_ylim(0, max(60, pct.values.max() + 8))
    ax.legend(fontsize=6, frameon=False)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.set_title("Composition (fold change above bars)", fontsize=8)
    p2 = os.path.join(FIG, "Fig_QC_orig_tsne_overlay_composition_preview.png")
    fig.savefig(p2, dpi=400, bbox_inches="tight", pad_inches=0.03); plt.close(fig); print("saved", p2)

    # 混合度（t-SNE 空间 30NN 的异样本邻居比例）
    nn = NearestNeighbors(n_neighbors=31).fit(U)
    _, idx = nn.kneighbors(U)
    other = np.array([(s[idx[i, 1:]] != s[i]).mean() for i in range(len(s))])
    print("\n=== 局部混合度（30NN 中异样本占比，0.5 = 完全混合）===")
    for g in ORD:
        m = ct == g
        print("  %-4s n=%5d | blank 细胞 %.2f | treated 细胞 %.2f" % (g, m.sum(), other[m & (s=="blank")].mean(), other[m & (s=="treated")].mean()))
    print("  总体 %.2f" % other.mean())
    pd.DataFrame({"ct": ct, "sample": s, "other_sample_frac": other}).to_csv(os.path.join(OUT, "tsne_mixing_by_cell.csv"), index=False)

    # 技术差异
    print("\n=== 两样本技术指标 ===")
    for smp in ("blank", "treated"):
        d = A.obs[A.obs["sample"] == smp]
        print("  %-8s n=%5d | UMI 中位 %6.0f | n_genes 中位 %5.0f | pct_mt 中位 %5.2f" % (
            smp, len(d), d.total_counts.median(), d.n_genes_by_counts.median(), d.pct_counts_mt.median()))

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\origQC_tsne.log", "w", encoding="utf8", buffering=1)
    main()
