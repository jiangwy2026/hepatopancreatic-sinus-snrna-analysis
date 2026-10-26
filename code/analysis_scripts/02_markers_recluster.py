import sys, io, os
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
from sklearn.metrics import adjusted_rand_score as ari
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sc.settings.verbosity = 1
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
base = root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result"

A = sc.read_h5ad(out + r"\adata_raw.h5ad")
cl = pd.read_csv(base + r"\2.Cluster\2.cluster\Cells.cluster.list.xls", sep="\t")
cl["Cells"] = cl["Cells"].astype(str)
A = A[A.obs_names.isin(set(cl["Cells"]))].copy()
mp = cl.set_index("Cells")
A.obs["vendor_cluster"] = mp.loc[A.obs_names, "Cluster"].astype(str).values
A.obs["group"] = mp.loc[A.obs_names, "Groups"].astype(str).values
print("cells:", A.shape[0], "| vendor clusters:", sorted(A.obs.vendor_cluster.unique(), key=int))
A.layers["counts"] = A.layers["counts"].copy()
A.X = A.layers["counts"].copy()
sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
Ann = A.copy()

# ---- marker genes (Wilcoxon) on vendor clusters ----
sc.tl.rank_genes_groups(Ann, "vendor_cluster", method="wilcoxon", pts=True)
res = []
for c in sorted(Ann.obs.vendor_cluster.unique(), key=int):
    d = sc.get.rank_genes_groups_df(Ann, group=c)
    d["cluster"] = c; res.append(d)
mk = pd.concat(res, ignore_index=True)
mk = mk.merge(A.var[["symbol"]].reset_index().rename(columns={"index":"names"}), on="names", how="left")
mk.to_csv(out + r"\markers_vendor_scanpy.csv", index=False)
print("\n=== 每群显著 marker 数 (padj<0.05 & |log2FC|>0.5) ===")
sig = mk[(mk.pvals_adj < 0.05) & (mk.logfoldchanges > 0.5)].groupby("cluster").size()
print(sig.to_string())
print("\n=== 每群 top8 marker（Wilcoxon, 按 score）===")
for c in sorted(Ann.obs.vendor_cluster.unique(), key=int):
    d = mk[mk.cluster == c].sort_values("scores", ascending=False).head(8)
    print("  C%-2s: %s" % (c, ", ".join(d.symbol.fillna(d.names).astype(str)[:8].tolist())))

# ---- 重聚类：多分辨率 ----
del Ann
sc.pp.highly_variable_genes(A, n_top_genes=2000, flavor="seurat_v3", layer="counts")
hv = A.var.highly_variable.values
B = A[:, hv].copy(); B.X = B.layers["counts"].copy()
sc.pp.normalize_total(B, target_sum=1e4); sc.pp.log1p(B)
sc.pp.scale(B, max_value=10)
sc.tl.pca(B, n_comps=50, svd_solver="arpack")
sc.pp.neighbors(B, n_neighbors=15, n_pcs=30)
rows = []
for r in [0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.5]:
    k = f"leiden_{r}"
    sc.tl.leiden(B, resolution=r, key_added=k, flavor="igraph", n_iterations=2, directed=False)
    a = ari(B.obs.vendor_cluster.astype(str).values, B.obs[k].astype(str).values)
    rows.append({"resolution": r, "n_clusters": B.obs[k].nunique(), "ARI_vs_vendor": round(a, 3)})
    print("  resol %.1f -> %d clusters, ARI vs vendor = %.3f" % (r, B.obs[k].nunique(), a))
pd.DataFrame(rows).to_csv(out + r"\leiden_resolution_summary.csv", index=False)
sc.tl.umap(B)
B.obs[["vendor_cluster","group","sample"]].to_csv(out + r"\umap_obs.csv")
pd.DataFrame(B.obsm["X_umap"], index=B.obs_names, columns=["UMAP1","UMAP2"]).to_csv(out + r"\umap_coords.csv")
B.write(out + r"\adata_hvg_processed.h5ad")
print("\nsaved markers_vendor_scanpy.csv, leiden_resolution_summary.csv, umap_coords.csv, adata_hvg_processed.h5ad")
