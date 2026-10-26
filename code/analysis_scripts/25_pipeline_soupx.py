import sys, io, os
import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
from sklearn.metrics import adjusted_rand_score as ari
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 250); sc.settings.verbosity = 1
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
A = sc.read_h5ad(out + r"\adata_soupx.h5ad")
A.X = A.layers["counts"].copy()
sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
sc.pp.highly_variable_genes(A, n_top_genes=2000, flavor="seurat_v3", layer="counts")
print("HVG:", int(A.var.highly_variable.sum()))
B = A[:, A.var.highly_variable].copy(); B.X = B.layers["counts"].copy()
sc.pp.normalize_total(B, target_sum=1e4); sc.pp.log1p(B); sc.pp.scale(B, max_value=10)
sc.tl.pca(B, n_comps=50, svd_solver="arpack")
sc.external.pp.harmony_integrate(B, "sample", basis="X_pca", adjusted_basis="X_pca_harmony")
sc.pp.neighbors(B, use_rep="X_pca_harmony", n_neighbors=15, n_pcs=30)
print("\n=== 校正后：多分辨率聚类（Harmony 整合）===")
rows = []
for r in [0.1,0.2,0.3,0.4,0.6,0.8,1.0]:
    k = "sx_%.1f" % r
    sc.tl.leiden(B, resolution=r, key_added=k, flavor="igraph", n_iterations=2, directed=False)
    a = ari(B.obs.vendor_cluster.astype(str).values, B.obs[k].astype(str).values)
    rows.append(dict(resolution=r, n_clusters=B.obs[k].nunique(), ARI_vs_原始分群=round(a,3)))
    print("  res %.1f -> %2d 群 | ARI vs 原始 = %.3f" % (r, B.obs[k].nunique(), a))
pd.DataFrame(rows).to_csv(out + r"\soupx_leiden_resolution.csv", index=False)
sc.tl.umap(B)
sc.settings.figdir = out + r"\figures_soupx"
sc.pl.umap(B, color=["vendor_cluster","sample"], show=False, save="_overview.png", ncols=2)
# 校正后 marker（按原始分群，便于前后对比）
sc.tl.rank_genes_groups(A, "vendor_cluster", method="wilcoxon", pts=True)
mk = pd.concat([sc.get.rank_genes_groups_df(A, group=g).assign(cluster=g) for g in A.obs.vendor_cluster.cat.categories], ignore_index=True)
mk["gene"] = mk["symbol"].fillna(mk["names"]) if "symbol" in mk else mk["names"]
mk = mk.merge(A.var.reset_index()[["index","sym"]].rename(columns={"index":"names"}), on="names", how="left")
mk["gene"] = mk["sym"].fillna(mk["names"]).astype(str)
mk.to_csv(out + r"\markers_soupx_vendor.csv", index=False)
sig = mk[(mk.pvals_adj<0.05) & (mk.logfoldchanges>0.5)]
print("\n=== 校正后：每群显著 marker 数（对比：校正前 111/219/258/201/370/58/469/219/5/6）===")
print(sig.groupby("cluster").size().to_string())
for g in sorted(A.obs.vendor_cluster.cat.categories, key=int):
    d = sig[sig.cluster==g].sort_values("scores", ascending=False).head(10)
    print("  C%-2s: %s" % (g, ", ".join(d.gene.astype(str).tolist())))
B.write(out + r"\adata_soupx_hvg.h5ad"); A.write(out + r"\adata_soupx_lognorm.h5ad")
print("\nsaved")
