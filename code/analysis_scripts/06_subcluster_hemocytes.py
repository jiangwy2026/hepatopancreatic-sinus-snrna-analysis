import sys, io, os
import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sc.settings.verbosity = 1
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
base = root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result"

A = sc.read_h5ad(out + r"\adata_raw.h5ad")
cl = pd.read_csv(base + r"\2.Cluster\2.cluster\Cells.cluster.list.xls", sep="\t"); cl["Cells"] = cl["Cells"].astype(str)
A = A[A.obs_names.isin(set(cl["Cells"]))].copy()
mp = cl.set_index("Cells"); A.obs["vendor_cluster"] = mp.loc[A.obs_names, "Cluster"].astype(str).values
A.obs["group"] = mp.loc[A.obs_names, "Groups"].astype(str).values
A.obs["sample"] = mp.loc[A.obs_names, "Samples"].astype(str).values

H = A[A.obs.vendor_cluster.isin(["0","1","2"])].copy()
print("血细胞 subset:", H.shape)
H.layers["counts"] = H.layers["counts"].copy()
H.X = H.layers["counts"].copy()
sc.pp.normalize_total(H, target_sum=1e4); sc.pp.log1p(H)
sc.pp.highly_variable_genes(H, n_top_genes=2000, flavor="seurat_v3", layer="counts")
print("HVG:", int(H.var.highly_variable.sum()))
Hs = H[:, H.var.highly_variable].copy(); Hs.X = Hs.layers["counts"].copy()
sc.pp.normalize_total(Hs, target_sum=1e4); sc.pp.log1p(Hs); sc.pp.scale(Hs, max_value=10)
sc.tl.pca(Hs, n_comps=50, svd_solver="arpack")
sc.pp.neighbors(Hs, n_neighbors=15, n_pcs=30)
from sklearn.metrics import silhouette_score
X = Hs.obsm["X_pca"][:, :30]
for r in [0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0]:
    k = f"leiden_{r}"
    sc.tl.leiden(Hs, resolution=r, key_added=k, flavor="igraph", n_iterations=2, directed=False)
    nl = Hs.obs[k].nunique()
    sizes = Hs.obs[k].value_counts()
    try: sil = silhouette_score(X, Hs.obs[k].astype(str).values, sample_size=3000, random_state=0)
    except Exception: sil = np.nan
    print("  res %.1f -> %2d subclusters (min size %3d, max %4d)  silhouette=%.3f" % (r, nl, sizes.min(), sizes.max(), sil))
sc.tl.umap(Hs)
Hs.write(out + r"\adata_hemocyte_sub.h5ad")
print("\nsaved adata_hemocyte_sub.h5ad")
