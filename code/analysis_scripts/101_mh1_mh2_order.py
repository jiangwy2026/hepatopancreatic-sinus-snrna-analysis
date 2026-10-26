# -*- coding: utf-8 -*-
"""
101_mh1_mh2_order.py - is the monocytic trajectory PH -> MH2 -> MH1, or is MH1 earlier?
Tests: (1) group position along the PH-rooted pseudotime, (2) CytoTRACE-style potential,
(3) when do MH2-up vs MH1-up genes peak along pseudotime (marker dynamics).
"""
import os, sys, json
import numpy as np, pandas as pd
import scanpy as sc, anndata as ad
import scipy.sparse as sp
from scipy.stats import mannwhitneyu

OUT = r"D:\桌面\single cell\analysis_output"
BASE = {"2": "PH", "1": "GH", "0": "MH1", "3": "MH2"}
ORD = ["PH", "GH", "MH1", "MH2"]

def main():
    A = ad.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))
    G = ad.read_h5ad(os.path.join(OUT, "adata_hemocyte_traj_geodesic.h5ad"))
    ct_all = A.obs["vendor_cluster"].astype(str).map(BASE)
    keep = ct_all.isin(ORD).values
    H = A[keep].copy()
    H.obs["ct"] = pd.Categorical(ct_all[keep].values, categories=ORD)
    assert list(H.obs_names) == list(G.obs_names), "cell order mismatch"
    H.obs["pt"] = np.asarray(G.obs["geodesic_pseudotime"], dtype=float)
    H.obs["pt_rank"] = pd.Series(H.obs["pt"]).rank(pct=True).values
    H.X = H.layers["counts"].copy() if "counts" in H.layers else H.X
    sc.pp.normalize_total(H, target_sum=1e4); sc.pp.log1p(H)
    print("cells:", H.shape, flush=True)

    # (1) where does each group sit on the pseudotime axis
    print("\n--- (1) pseudotime position (raw geodesic; rank) ---", flush=True)
    t1 = H.obs.groupby("ct", observed=True)["pt"].agg(["count", "median", "mean"])
    t1["rank_median"] = H.obs.groupby("ct", observed=True)["pt_rank"].median()
    print(t1.round(3).to_string(), flush=True)

    # (2) CytoTRACE-style potential (smoothed n_genes)
    ng = np.asarray(H.obs["n_genes_by_counts"], dtype=float)
    H.raw = None
    W = None
    try:
        A2 = H.copy(); sc.pp.highly_variable_genes(A2, n_top_genes=2000); A2 = A2[:, A2.var.highly_variable]
        sc.pp.scale(A2, max_value=10); sc.tl.pca(A2, n_comps=30); sc.pp.neighbors(A2, n_neighbors=30)
        W = A2.obsp["connectivities"].tocsr()
    except Exception as e:
        print("neighbour rebuild failed:", repr(e)[:80], flush=True)
    if W is not None:
        rs = np.asarray(W.sum(1)).ravel(); rs[rs == 0] = 1
        gcs = np.asarray(W.dot(ng) / rs)
    else:
        gcs = ng
    print("\n--- (2) CytoTRACE-style score by group (higher = more stem-like) ---", flush=True)
    print(pd.DataFrame({"ct": H.obs["ct"].astype(str).values, "gcs": gcs}).groupby("ct")["gcs"].mean().round(1).to_string(), flush=True)

    # (3) marker dynamics: genes up in MH2 vs up in MH1, where do they peak?
    H.obs["ct"] = pd.Categorical(H.obs["ct"].astype(str), categories=ORD)
    sc.tl.rank_genes_groups(H, "ct", groups=["MH2"], reference="MH1", method="wilcoxon", n_genes=400)
    up_mh2 = [g for g in H.uns["rank_genes_groups"]["names"]["MH2"][:150]]
    sc.tl.rank_genes_groups(H, "ct", groups=["MH1"], reference="MH2", method="wilcoxon", n_genes=400)
    up_mh1 = [g for g in H.uns["rank_genes_groups"]["names"]["MH1"][:150]]
    def dense(names):
        idx = [H.var_names.get_loc(g) for g in names if g in H.var_names]
        X = H.X[:, idx]
        return np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X), [H.var_names[i] for i in idx]
    def peak_dyn(mat):
        q = np.quantile(H.obs["pt_rank"].values, np.linspace(0, 1, 21))
        b = np.clip(np.digitize(H.obs["pt_rank"].values, q) - 1, 0, 19)
        out = []
        for j in range(mat.shape[1]):
            prof = np.array([mat[b == i, j].mean() for i in range(20)])
            if prof.max() <= 0:
                continue
            k = int(np.argmax(prof))
            com = float((prof * np.arange(20)).sum() / prof.sum())
            out.append((k, com))
        return out
    print("\n--- (3) marker dynamics along pseudotime (bin 0 = earliest) ---", flush=True)
    res = {}
    for label, genes in [("MH2-up", up_mh2), ("MH1-up", up_mh1)]:
        m, used = dense(genes)
        dyn = peak_dyn(m)
        ks = np.array([d[0] for d in dyn]); com = np.array([d[1] for d in dyn])
        res[label] = {"genes": used, "peak_bins": ks.tolist(), "com": com.tolist()}
        print("   %-7s n=%d | median peak bin %.1f | median centre-of-mass %.2f" % (
            label, len(ks), np.median(ks), np.median(com)), flush=True)
    if "MH2-up" in res and "MH1-up" in res:
        u = mannwhitneyu(res["MH2-up"]["com"], res["MH1-up"]["com"], alternative="less")
        print("   Mann-Whitney (MH2-up earlier than MH1-up): U=%.0f  p=%.3g" % (u.statistic, u.pvalue), flush=True)
        u2 = mannwhitneyu(res["MH2-up"]["peak_bins"], res["MH1-up"]["peak_bins"], alternative="less")
        print("   Mann-Whitney on peak bins:                U=%.0f  p=%.3g" % (u2.statistic, u2.pvalue), flush=True)
    pd.DataFrame({"ct": H.obs["ct"].astype(str).values, "pt": H.obs["pt"].values}).to_csv(os.path.join(OUT, "traj_group_pseudotime.csv"), index=False)
    json.dump({"group_pseudotime": t1.to_dict(), "gcs": pd.DataFrame({"ct": H.obs["ct"].astype(str).values, "gcs": gcs}).groupby("ct")["gcs"].mean().to_dict(),
               "mh2_up_first20": up_mh2[:20], "mh1_up_first20": up_mh1[:20]},
              open(os.path.join(OUT, "MH1_MH2_ORDER_REPORT.json"), "w"), indent=1)

if __name__ == "__main__":
    sys.stdout = open(r"D:\桌面\single cell\logs\mh1_mh2_order.log", "w", encoding="utf8", buffering=1)
    main()
