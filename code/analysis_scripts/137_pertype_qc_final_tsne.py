# -*- coding: utf-8 -*-
"""137_pertype_qc_final_tsne.py - 修正版：阈值由"干净细胞"按类型估计，再应用到全部细胞"""
import io, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
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

def main():
    A = ad.read_h5ad(os.path.join(OUT, "adata_raw.h5ad"))
    ct = pd.read_csv(os.path.join(OUT, "all_19795_celltype.csv"), index_col=0)["ct"].astype(str)
    A.obs["ct"] = pd.Categorical(ct.reindex(list(map(str, A.obs_names))).values, categories=ORD)
    sx = ad.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"), backed="r")
    keptset = set(map(str, sx.obs_names)); del sx
    A.obs["clean"] = [str(x) in keptset for x in A.obs_names]
    A.X = A.layers["counts"].copy()

    keep = np.zeros(A.n_obs, bool); rows = []
    for g in ORD:
        m = (A.obs["ct"] == g).values
        clean = m & A.obs["clean"].values
        mt_c = A.obs.pct_counts_mt.values[clean]
        med = np.median(mt_c); mad = np.median(np.abs(mt_c - med)) * 1.4826
        cut = med + 3 * mad
        k = m & (A.obs.pct_counts_mt.values < cut) & (A.obs.n_genes_by_counts.values >= 200)
        keep |= k
        rows.append({"ct": g, "n_all": int(m.sum()), "mt_med_clean": round(med, 2), "cut": round(cut, 2),
                     "kept": int(k.sum()), "added": int((k & ~A.obs["clean"].values).sum())})
    T = pd.DataFrame(rows); print(T.to_string(index=False))
    print("合计保留 %d（厂商 %d，新增 %d）" % (keep.sum(), int(A.obs["clean"].sum()), int((keep & ~A.obs["clean"].values).sum())))
    T.to_csv(os.path.join(OUT, "pertype_thresholds_final.csv"), index=False)
    B = A[keep].copy(); B.obs["added"] = ~B.obs["clean"].values

    import scanpy.external as sce
    sc.pp.normalize_total(B, target_sum=1e4); sc.pp.log1p(B)
    sc.pp.highly_variable_genes(B, n_top_genes=2000)
    Bh = B[:, B.var.highly_variable].copy()
    sc.pp.scale(Bh, max_value=10); sc.tl.pca(Bh, n_comps=30, svd_solver="arpack")
    sce.pp.harmony_integrate(Bh, "sample")
    sc.pp.neighbors(Bh, n_neighbors=15, use_rep="X_pca_harmony")
    sc.tl.tsne(Bh, n_pcs=30, perplexity=30, random_state=0, n_jobs=8)
    B.obsm["X_tsne"] = Bh.obsm["X_tsne"]
    B.write_h5ad(os.path.join(OUT, "adata_pertypeQC_final_tsne.h5ad"))
    U = np.asarray(B.obsm["X_tsne"]); s = B.obs["sample"].astype(str).values
    c = B.obs["ct"].astype(str).values

    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.9))
    for ax, smp, title in zip(axes, ("blank", "treated"), ("Control (0 h)", "Hemospasia (6 h)")):
        m = s == smp
        for g in ORD:
            mm = m & (c == g)
            if mm.sum():
                ax.scatter(U[mm, 0], U[mm, 1], s=1.2, c=COL[g], linewidths=0, rasterized=True,
                           label="%s (%d)" % (g, int(mm.sum())))
        ax.set_title("%s   n = %s" % (title, format(int(m.sum()), ",")), fontsize=8)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("t-SNE 1"); ax.set_ylabel("t-SNE 2"); ax.set_box_aspect(1)
        for sp_ in ax.spines.values(): sp_.set_color("#999999")
    axes[1].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, markerscale=4,
                   fontsize=6.2, handletextpad=0.3, labelspacing=0.3)
    fig.suptitle("Per-cell-type QC (thresholds from clean cells: pct_mt < median + 3xMAD) - %s cells" % format(len(s), ","),
                 fontsize=8, y=1.0)
    p = os.path.join(FIG, "Fig_QC_pertype_tsne_blank_treated_preview.png")
    fig.savefig(p, dpi=400, bbox_inches="tight", pad_inches=0.03); plt.close(fig); print("saved", p)

    # 组成对比
    t = pd.crosstab(B.obs["ct"], B.obs["sample"])
    pc = (t / t.sum(axis=0) * 100).reindex(ORD)
    print("\n=== 修正方案下的组成 ===")
    print(pc.round(2).to_string())
    print("倍数 (treated/blank):", (pc["treated"] / pc["blank"]).round(2).to_dict())
    pc.to_csv(os.path.join(OUT, "composition_pertypeQC_final.csv"))

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\pertype_tsne_final.log", "w", encoding="utf8", buffering=1)
    main()
