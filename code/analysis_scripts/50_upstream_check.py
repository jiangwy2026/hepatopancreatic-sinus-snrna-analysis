import sys, io, os
import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
from scipy.stats import rankdata, spearmanr
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 250); sc.settings.verbosity = 0
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
A = sc.read_h5ad(out + r"\adata_soupx_hvg.h5ad")   # 全细胞 + Harmony 坐标
F = sc.read_h5ad(out + r"\adata_soupx.h5ad")
# ---- 1) CytoTRACE 式"分化潜能"打分（在全部 9658 细胞上）----
F.X = F.layers["counts"].copy(); sc.pp.normalize_total(F, target_sum=1e4); sc.pp.log1p(F)
gcs = np.asarray((F.layers["counts"] > 0).sum(1)).ravel().astype(float)
X = F.X
rng = np.random.default_rng(0); si = rng.choice(F.n_obs, size=min(6000, F.n_obs), replace=False)
rg = rankdata(gcs[si]); rg = (rg-rg.mean())/(rg.std()+1e-12)
var = np.asarray(X[si].multiply(X[si]).mean(0)).ravel() - np.asarray(X[si].mean(0)).ravel()**2
hv = np.argsort(-var)[:6000]
Xs = X[si][:, hv].toarray().astype(np.float32)
rx = rankdata(Xs, axis=0); rx = (rx-rx.mean(0))/(rx.std(0)+1e-12)
rho_g = (rx*rg[:,None]).mean(0)
top = np.argsort(-rho_g)[:200]
F.obs["cytotrace_like"] = np.asarray(X[:, hv[top]].mean(1)).ravel()
print("=== CytoTRACE 式分化潜能（越高=越原始）按群平均 ===")
g = F.obs.groupby(F.obs.vendor_cluster.astype(str)).agg(n=("cytotrace_like","size"), 潜能为=("cytotrace_like","mean"),
        中位=("cytotrace_like","median"), 每细胞基因数=("n_genes_by_counts","median"))
print(g.sort_values("潜能为", ascending=False).round(3).to_string())
print("\n=== 潜能 vs 细胞周期（G2M）相关（血细胞群内）===")
H = F[F.obs.vendor_cluster.astype(str).isin(["0","1","2"])]
print("  Spearman(潜能为, n_genes) = %.3f" % spearmanr(H.obs["cytotrace_like"], H.obs["n_genes_by_counts"])[0])
# ---- 2) 全局 PAGA（校正后）看谁连到 C2 ----
A2 = A.copy(); A2.obs["grp"] = "C" + A2.obs.vendor_cluster.astype(str)
sc.pp.neighbors(A2, use_rep="X_pca_harmony", n_neighbors=15, n_pcs=30)
sc.tl.paga(A2, groups="grp")
ct = pd.DataFrame(A2.uns["paga"]["connectivities"].todense(), index=A2.obs["grp"].cat.categories, columns=A2.obs["grp"].cat.categories)
print("\n=== 全局 PAGA（校正后）— C2 与其他群的连接强度 ===")
print(ct.loc["C2"].drop("C2").sort_values(ascending=False).round(3).to_string())
# ---- 3) C2 内部是否有更原始的亚结构 ----
H2 = F[F.obs.vendor_cluster.astype(str) == "2"].copy()
sc.pp.highly_variable_genes(H2, n_top_genes=2000, flavor="seurat_v3", layer="counts")
Hb = H2[:, H2.var.highly_variable].copy(); Hb.X = Hb.layers["counts"].copy()
sc.pp.normalize_total(Hb, target_sum=1e4); sc.pp.log1p(Hb); sc.pp.scale(Hb, max_value=10)
sc.tl.pca(Hb, n_comps=30, svd_solver="arpack")
sc.pp.neighbors(Hb, n_neighbors=15, n_pcs=20)
sc.tl.leiden(Hb, resolution=0.4, key_added="c2sub", flavor="igraph", n_iterations=2, directed=False)
print("\n=== C2 内部亚结构（res 0.4）===")
print(Hb.obs["c2sub"].value_counts().sort_index().to_string())
F.obs["cyto"] = F.obs["cytotrace_like"]
sub = F.loc[Hb.obs_names].copy(); sub.obs["c2sub"] = Hb.obs["c2sub"].values
print(sub.obs.groupby("c2sub").agg(n=("cyto","size"), 潜能=("cyto","mean"), G2M=("G2M_score","mean" if "G2M_score" in sub.obs else "size")).round(3).to_string() if "G2M_score" in sub.obs else sub.obs.groupby("c2sub").agg(n=("cyto","size"), 潜能=("cyto","mean")).round(3).to_string())
F.write(out + r"\adata_soupx_with_potency.h5ad")
print("\nsaved adata_soupx_with_potency.h5ad")
