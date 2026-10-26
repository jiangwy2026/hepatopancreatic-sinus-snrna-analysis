import sys, io, os
import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, scanpy as sc
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 260); sc.settings.verbosity = 0
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
R = sc.read_h5ad(out + r"\geo_GSE306674_raw.h5ad")
print("原始:", R.shape)
R = R[R.obs.n_genes_by_counts >= 200].copy()
sc.pp.filter_genes(R, min_cells=3)
print("过滤后:", R.shape)
R.X = R.layers["counts"].copy()
sc.pp.normalize_total(R, target_sum=1e4); sc.pp.log1p(R)
sc.pp.highly_variable_genes(R, n_top_genes=2000, flavor="seurat_v3", layer="counts")
Rh = R[:, R.var.highly_variable].copy(); Rh.X = Rh.layers["counts"].copy()
sc.pp.normalize_total(Rh, target_sum=1e4); sc.pp.log1p(Rh); sc.pp.scale(Rh, max_value=10)
sc.tl.pca(Rh, n_comps=50, svd_solver="arpack")
sc.external.pp.harmony_integrate(Rh, "sample", basis="X_pca", adjusted_basis="X_pca_harmony")
sc.pp.neighbors(Rh, use_rep="X_pca_harmony", n_neighbors=15, n_pcs=30)
sc.tl.leiden(Rh, resolution=0.5, key_added="cl", flavor="igraph", n_iterations=2, directed=False)
print("\n=== 参考数据聚类（res 0.5）===")
print(Rh.obs["cl"].value_counts().sort_index().to_string())
sc.tl.rank_genes_groups(R, "cl", groups=[g for g in Rh.obs["cl"].cat.categories], method="wilcoxon", pts=True) if False else None
Rh = Rh.raw if False else Rh
sc.tl.rank_genes_groups(Rh, "cl", method="wilcoxon", pts=True)
mk = pd.concat([sc.get.rank_genes_groups_df(Rh, group=g).assign(cl=g) for g in Rh.obs["cl"].cat.categories], ignore_index=True)
mk.to_csv(out + r"\geo_ref_cluster_markers.csv", index=False)
print("\n=== 各参考群 top10 marker ===")
for g in Rh.obs["cl"].cat.categories:
    d = mk[mk.cl == g].sort_values("scores", ascending=False).head(10)
    n = int((Rh.obs["cl"] == g).sum())
    print("  R%-3s (n=%5d): %s" % (g, n, ", ".join(d["names"].astype(str).tolist())))
Rh.write(out + r"\geo_GSE306674_processed.h5ad")
print("\nsaved")
