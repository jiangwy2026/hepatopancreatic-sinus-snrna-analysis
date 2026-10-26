import sys, io, os, traceback
import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, scanpy as sc
from sklearn.metrics import silhouette_score
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 250); sc.settings.verbosity = 1
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
base = root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result"

A = sc.read_h5ad(out + r"\adata_raw.h5ad")
cl = pd.read_csv(base + r"\2.Cluster\2.cluster\Cells.cluster.list.xls", sep="\t"); cl["Cells"]=cl["Cells"].astype(str)
A = A[A.obs_names.isin(set(cl["Cells"]))].copy()
mp = cl.set_index("Cells"); A.obs["vendor_cluster"]=mp.loc[A.obs_names,"Cluster"].astype(str).values
A.obs["group"]=mp.loc[A.obs_names,"Groups"].astype(str).values
A.obs["sample"]=mp.loc[A.obs_names,"Samples"].astype(str).values
H = A[A.obs.vendor_cluster.isin(["0","1","2"])].copy()
H.layers["counts"]=H.layers["counts"].copy(); H.X=H.layers["counts"].copy()
sc.pp.normalize_total(H, target_sum=1e4); sc.pp.log1p(H)
sc.pp.highly_variable_genes(H, n_top_genes=2000, flavor="seurat_v3", layer="counts")
H = H[:, H.var.highly_variable].copy(); H.X = H.layers["counts"].copy()
sc.pp.normalize_total(H, target_sum=1e4); sc.pp.log1p(H); sc.pp.scale(H, max_value=10)
sc.tl.pca(H, n_comps=50, svd_solver="arpack")

# Harmony 整合（按样本）
try:
    sc.external.pp.harmony_integrate(H, "sample", basis="X_pca", adjusted_basis="X_pca_harmony")
except Exception:
    try:
        sc.pp.harmony_integrate(H, "sample", basis="X_pca", adjusted_basis="X_pca_harmony")
    except Exception:
        traceback.print_exc(); raise
print("harmony done:", H.obsm["X_pca_harmony"].shape)

print("\n=== 整合后：各分辨率结果 ===")
sc.pp.neighbors(H, n_neighbors=15, n_pcs=30, use_rep="X_pca_harmony", key_added="hn")
X = H.obsm["X_pca_harmony"][:, :30]
for r in [0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0]:
    k=f"sub_{r}"
    sc.tl.leiden(H, resolution=r, key_added=k, neighbors_key="hn", flavor="igraph", n_iterations=2, directed=False)
    s=H.obs[k].value_counts()
    try: sil=silhouette_score(X, H.obs[k].astype(str).values, sample_size=3000, random_state=0)
    except Exception: sil=float("nan")
    print("  res %.1f -> %2d 群 (min %3d, max %4d) silhouette=%.3f" % (r, len(s), s.min(), s.max(), sil))
H.write(out + r"\adata_hemocyte_harmony.h5ad")
print("\nsaved adata_hemocyte_harmony.h5ad")
