import sys, io, os, re
import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, scanpy as sc, anndata as ad, scipy.sparse as sp
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 250); sc.settings.verbosity = 0
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
base = root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result"
# 我的数据（未校正，与 eLife 预处理匹配）
A = sc.read_h5ad(out + r"\adata_raw.h5ad")
cl = pd.read_csv(base + r"\2.Cluster\2.cluster\Cells.cluster.list.xls", sep="\t"); cl["Cells"]=cl["Cells"].astype(str)
A = A[A.obs_names.isin(set(cl["Cells"]))].copy()
A.obs["my_cluster"] = cl.set_index("Cells").reindex(A.obs_names)["Cluster"].astype(str).values
A.obs["dataset"] = "Mine"
A.obs["sub"] = "Mine_C" + A.obs["my_cluster"].astype(str)
E = sc.read_h5ad(out + r"\elife_adata.h5ad")
E.obs["dataset"] = "eLife"
E.obs["sub"] = "eLife_C" + E.obs["elife_cluster"].astype(str)
shared = sorted(set(A.var_names) & set(E.var_names))
print("共享基因:", len(shared), "| 我的 %d | eLife %d" % (A.n_vars, E.n_vars))
A = A[:, shared].copy(); E = E[:, shared].copy()
C = ad.concat([A, E], join="inner", index_unique=None, label=None, keys=None, merge="same")
C.obs["dataset"] = C.obs["dataset"].astype(str)
C.layers["counts"] = C.X.copy()
print("合并:", C.shape)
print(C.obs.dataset.value_counts().to_string())
C.X = C.layers["counts"].copy()
sc.pp.normalize_total(C, target_sum=1e4); sc.pp.log1p(C)
sc.pp.highly_variable_genes(C, n_top_genes=2000, flavor="seurat_v3", layer="counts")
H = C[:, C.var.highly_variable].copy(); H.X = H.layers["counts"].copy()
sc.pp.normalize_total(H, target_sum=1e4); sc.pp.log1p(H); sc.pp.scale(H, max_value=10)
sc.tl.pca(H, n_comps=50, svd_solver="arpack")
sc.external.pp.harmony_integrate(H, "dataset", basis="X_pca", adjusted_basis="X_pca_harmony")
sc.pp.neighbors(H, use_rep="X_pca_harmony", n_neighbors=15, n_pcs=30)
sc.tl.umap(H)
sc.tl.leiden(H, resolution=0.5, key_added="cl_joint", flavor="igraph", n_iterations=2, directed=False)
H.write(out + r"\joint_mine_elife.h5ad")
print("联合聚类:", H.obs.cl_joint.nunique(), "群")
# ---- 指标 1: 混入度 ----
import scipy.sparse as sp2
K = H.obsp["connectivities"].tocoo()
ds = (H.obs["dataset"] == "eLife").values.astype(int)
same = ds[K.row] == ds[K.col]
cross = (~same).sum() / len(K.row)
print("\n=== 整体混入度 ===")
print("  跨数据集连边比例: %.1f%% (50%% 为完全随机混合)" % (cross*100))
# 每个联合群的数据集组成
t = pd.crosstab(H.obs.cl_joint, H.obs.dataset)
t["我的占比%"] = (t["Mine"]/t.sum(axis=1)*100).round(1)
print("\n=== 各联合群的数据集组成 ===")
print(t.to_string())
# ---- 指标 2: 我的细胞在 UMAP 上的邻居（谁是 eLife）----
mine = (H.obs.dataset == "Mine").values
el   = (H.obs.dataset == "eLife").values
indptr = H.obsp["connectivities"].tocsr()
rows = H.obsp["connectivities"].tocsr()
mine_idx = np.where(mine)[0]
res = []
elc = H.obs["elife_cluster"].astype(str).values
for i in mine_idx:
    nb = rows.indices[rows.indptr[i]:rows.indptr[i+1]]
    if len(nb) == 0: continue
    enb = nb[el[nb]]
    res.append({"my_cluster": H.obs["my_cluster"].values[i],
                "n_nb": len(nb), "n_elife_nb": len(enb),
                "frac_elife": len(enb)/len(nb),
                "top_elife_cluster": pd.Series(elc[enb]).mode().iloc[0] if len(enb) else None})
R = pd.DataFrame(res)
print("\n=== 我的每个群：邻居中 eLife 细胞的比例与最常邻接的 eLife 群 ===")
g = R.groupby("my_cluster").agg(n=("frac_elife","size"), 邻居含eLife比例=("frac_elife","mean"),
                                最常邻接eLife群=("top_elife_cluster", lambda s: s.mode().iloc[0] if len(s.mode()) else "NA"))
print(g.round(3).to_string())
R.to_csv(out + r"\joint_neighbor_composition.csv", index=False)
