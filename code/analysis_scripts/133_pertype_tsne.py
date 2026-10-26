# -*- coding: utf-8 -*-
"""
133_pertype_tsne.py - 按类型阈值（median+3MAD）保留细胞后重做降维，出 control / treated 两张 t-SNE
"""
import io, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
import scipy.sparse as sp
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
sc.settings.verbosity = 1
OUT = r"D:\桌面\single cell\analysis_output"; FIG = os.path.join(OUT, "figures_pub")
ORD = ["PH", "GH", "MH1", "MH2", "Epi", "F", "R", "B", "NEU"]
COL = {"PH": "#B2182B", "GH": "#E08214", "MH1": "#FDB863", "MH2": "#8C510A", "Epi": "#7B3294",
       "F": "#006D2C", "R": "#41AB5D", "B": "#78C679", "NEU": "#2166AC"}
RULE = float(os.environ.get("K_MAD", "3"))

def main():
    A = ad.read_h5ad(os.path.join(OUT, "adata_raw.h5ad"))
    ct = pd.read_csv(os.path.join(OUT, "all_19795_celltype.csv"), index_col=0)["ct"].astype(str)
    A.obs["ct"] = pd.Categorical(ct.reindex(list(map(str, A.obs_names))).values, categories=ORD)
    _sx = ad.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"), backed="r")
    _keptset = set(map(str, _sx.obs_names))
    A.obs["orig_kept"] = [str(x) in _keptset for x in A.obs_names]
    del _sx
    A.X = A.layers["counts"].copy()
    # ---- 按类型阈值
    keep = np.zeros(A.n_obs, bool)
    info = []
    for g in ORD:
        m = (A.obs["ct"] == g).values
        mt = A.obs.pct_counts_mt.values[m]; ng = A.obs.n_genes_by_counts.values[m]
        med = np.median(mt); mad = np.median(np.abs(mt - med)) * 1.4826
        cut = med + RULE * mad
        k = (mt < cut) & (ng >= 200)
        keep[m] = k
        info.append({"ct": g, "n_all": int(m.sum()), "mt_med": round(med, 1), "cut": round(cut, 1),
                     "kept": int(k.sum()), "added_vs_vendor": int((k & ~A.obs.orig_kept.values[m]).sum())})
    I = pd.DataFrame(info); print(I.to_string(index=False))
    print("保留合计 %d (原厂商 %d) | 新增 %d" % (keep.sum(), int(A.obs.orig_kept.sum()), int((keep & ~A.obs.orig_kept.values).sum())))
    B = A[keep].copy()
    I.to_csv(os.path.join(OUT, "pertype_thresholds_used.csv"), index=False)
    B.obs["added"] = ~B.obs["orig_kept"].values

    # ---- 降维（Harmony 按样本校正）
    import scanpy.external as sce
    sc.pp.normalize_total(B, target_sum=1e4); sc.pp.log1p(B)
    sc.pp.highly_variable_genes(B, n_top_genes=2000)
    Bh = B[:, B.var.highly_variable].copy()
    sc.pp.scale(Bh, max_value=10); sc.tl.pca(Bh, n_comps=30, svd_solver="arpack")
    sce.pp.harmony_integrate(Bh, "sample")
    sc.pp.neighbors(Bh, n_neighbors=15, use_rep="X_pca_harmony")
    sc.tl.tsne(Bh, n_pcs=30, perplexity=30, random_state=0, n_jobs=8)
    B.obsm["X_tsne"] = Bh.obsm["X_tsne"]
    B.write_h5ad(os.path.join(OUT, "adata_pertypeQC_harmony_tsne.h5ad"))
    U = np.asarray(B.obsm["X_tsne"]); s = B.obs["sample"].astype(str).values; c = B.obs["ct"].astype(str).values

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.9))
    for ax, smp, title in zip(axes, ("blank", "treated"), ("Control (0 h)", "Hemospasia 6 h")):
        m = s == smp
        for g in ORD:
            mm = m & (c == g)
            if mm.sum() == 0: continue
            ax.scatter(U[mm, 0], U[mm, 1], s=1.2, c=COL[g], linewidths=0, rasterized=True,
                       label="%s (%d)" % (g, int(mm.sum())))
        ax.set_title("%s   n=%d" % (title, int(m.sum())), fontsize=8)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_box_aspect(1)
        ax.set_xlabel("t-SNE 1"); ax.set_ylabel("t-SNE 2")
        for sp_ in ax.spines.values(): sp_.set_color("#999999")
    axes[1].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, markerscale=4,
                   fontsize=6.2, handletextpad=0.3, labelspacing=0.3)
    fig.suptitle("Per-cell-type QC (pct_mt < median + %gxMAD, n_genes >= 200); %d cells shown" % (RULE, B.n_obs),
                 fontsize=8, y=1.0)
    p = os.path.join(FIG, "Fig_QC_pertype_tsne_blank_treated_preview.png")
    fig.savefig(p, dpi=400, bbox_inches="tight", pad_inches=0.03); plt.close(fig)
    print("saved", p)

    # ---- 诊断：新加入的细胞落在哪里
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.9))
    for ax, smp in zip(axes[:2], ("blank", "treated")):
        m = s == smp
        ax.scatter(U[m & ~B.obs["added"].values, 0], U[m & ~B.obs["added"].values, 1], s=1.2,
                   c="#D9D9D9", linewidths=0, rasterized=True, label="原厂商保留")
        ax.scatter(U[m & B.obs["added"].values, 0], U[m & B.obs["added"].values, 1], s=1.4,
                   c="#2166AC", linewidths=0, rasterized=True, label="按类型阈值新纳入")
        ax.set_title("%s: 新纳入 %d 个" % (smp, int((m & B.obs["added"].values).sum())), fontsize=8)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_box_aspect(1); ax.legend(fontsize=6, frameon=False, markerscale=4)
        for sp_ in ax.spines.values(): sp_.set_color("#999999")
    ax = axes[2]
    ax.scatter(U[~B.obs["added"].values, 0], U[~B.obs["added"].values, 1], s=1.2, c="#D9D9D9", linewidths=0, rasterized=True)
    mm = B.obs["added"].values & (B.obs.pct_counts_mt.values > 20)
    ax.scatter(U[mm, 0], U[mm, 1], s=1.4, c="#B2182B", linewidths=0, rasterized=True, label="新纳入且 pct_mt>20%")
    ax.set_title("新纳入细胞中高线粒体者 (n=%d)" % int(mm.sum()), fontsize=8)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_box_aspect(1); ax.legend(fontsize=6, frameon=False, markerscale=4)
    p2 = os.path.join(FIG, "Fig_QC_pertype_added_cells_preview.png")
    fig.savefig(p2, dpi=400, bbox_inches="tight", pad_inches=0.03); plt.close(fig)
    print("saved", p2)

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\pertype_tsne.log", "w", encoding="utf8", buffering=1)
    main()
