# -*- coding: utf-8 -*-
"""
93_traj_rebuild.py - rebuild the hemocyte trajectory with a proper diffusion map
(the earlier version used Harmony-corrected PCs, which produced no spectral gap) and
test whether a branching (Y-shaped) structure like the reference figure is recoverable.
"""
import os, sys, json
import numpy as np, pandas as pd
import scanpy as sc, anndata as ad
import scipy.sparse as sp

OUT = r"D:\桌面\single cell\analysis_output"
BASE = {"2": "PH", "1": "GH", "0": "MH1", "3": "MH2"}
BLOOD = ["PH", "GH", "MH1", "MH2"]

def prep():
    A = ad.read_h5ad(os.path.join(OUT, "adata_soupx_h5ad".replace("_h5ad", ".h5ad") if False else "adata_soupx.h5ad"))
    ct = A.obs["vendor_cluster"].astype(str).map(BASE)
    keep = ct.isin(BLOOD).values
    H = A[keep].copy()
    H.obs["ct"] = pd.Categorical(ct[keep].values, categories=BLOOD)
    H.X = H.layers["counts"].copy() if "counts" in H.layers else H.X
    sc.pp.normalize_total(H, target_sum=1e4)
    sc.pp.log1p(H)
    H.raw = H
    sc.pp.highly_variable_genes(H, n_top_genes=2000)
    H = H[:, H.var.highly_variable].copy()
    sc.pp.scale(H, max_value=10)
    print("subset:", H.shape, "| groups:", H.obs["ct"].value_counts().to_dict(), flush=True)
    return H

def run(H, tag, n_pcs=30, n_neighbors=30, harmony=False):
    X = H.copy()
    sc.tl.pca(X, n_comps=n_pcs, svd_solver="arpack")
    pca_key = "X_pca"
    if harmony:
        try:
            sce.pp.harmony_integrate(X, "sample")
            pca_key = "X_pca_harmony"
        except Exception as e:
            print(tag, "harmony failed:", repr(e)[:80]); return None
    sc.pp.neighbors(X, n_neighbors=n_neighbors, use_rep=pca_key)
    sc.tl.diffmap(X, n_comps=15)
    ev = np.asarray(X.uns["diffmap_evals"])
    gaps = ev[:-1] - ev[1:]
    dc = np.asarray(X.obsm["X_diffmap"])[:, 1:4]
    ct = X.obs["ct"].astype(str).values
    # group separation in diffusion space (silhouette on DC2/DC3, which usually carry the branches)
    from sklearn.metrics import silhouette_score
    try:
        sil = silhouette_score(dc[:, :2], ct)
    except Exception:
        sil = float("nan")
    means = {g: np.round(dc[ct == g].mean(axis=0), 4).tolist() for g in BLOOD}
    print("[%s] n_pcs=%d nn=%d harmony=%s" % (tag, n_pcs, n_neighbors, harmony), flush=True)
    print("   evals:", np.round(ev, 4).tolist(), flush=True)
    print("   largest gaps:", np.round(np.sort(gaps)[-3:], 4).tolist(), flush=True)
    print("   mean DC2/DC3 by group:", means, flush=True)
    print("   silhouette(4 groups in DC2/DC3): %.3f" % sil, flush=True)
    # correlation of DC2 with technical covariates
    for cov in ("total_counts", "pct_counts_mt"):
        if cov in X.obs:
            r = np.corrcoef(dc[:, 0], np.asarray(X.obs[cov], dtype=float))[0, 1]
            print("   corr(DC2, %s) = %+.3f" % (cov, r), flush=True)
    r_s = np.corrcoef(dc[:, 0], (X.obs["sample"] == "treated").astype(float))[0, 1]
    print("   corr(DC2, treated) = %+.3f" % r_s, flush=True)
    X.uns["traj_tag"] = tag
    return X, sil, float(np.sort(gaps)[-1])

if __name__ == "__main__":
    import scanpy.external as sce
    sys.stdout = open(r"D:\桌面\single cell\logs\traj_rebuild.log", "w", encoding="utf8", buffering=1)
    H = prep()
    best = None
    for tag, kw in [("plain-nn30", dict(n_neighbors=30, harmony=False)),
                    ("plain-nn15", dict(n_neighbors=15, harmony=False)),
                    ("harmony-nn30", dict(n_neighbors=30, harmony=True))]:
        out = run(H, tag, **kw)
        if out is None: continue
        X, sil, gap = out
        score = sil
        if best is None or score > best[2]:
            best = (tag, X, score)
        X.write_h5ad(os.path.join(OUT, "adata_hemocyte_traj_%s.h5ad" % tag))
    print("\nbest by silhouette:", best[0], "%.3f" % best[2], flush=True)
