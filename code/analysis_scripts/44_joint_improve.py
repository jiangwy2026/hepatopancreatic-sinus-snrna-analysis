import sys, io, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 260); sc.settings.verbosity = 0
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
base = root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result"
H = sc.read_h5ad(out + r"\joint_mine_elife.h5ad")
mix_edges = lambda HH: (lambda K: 1 - (HH.obs["dataset"].values[K.row] == HH.obs["dataset"].values[K.col]).mean())(HH.obsp["connectivities"].tocoo())
print("=== 当前（2000 HVG, Harmony 默认）跨数据集连边 %.1f%% ===" % (mix_edges(H)*100))
# --- 加强版：3000 HVG、batch-aware HVG、回归深度、Harmony 20 轮 ---
A = H.copy()
A.X = A.layers["counts"].copy()
sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
sc.pp.highly_variable_genes(A, n_top_genes=3000, flavor="seurat_v3", layer="counts", batch_key="dataset")
B = A[:, A.var.highly_variable].copy(); B.X = B.layers["counts"].copy()
sc.pp.normalize_total(B, target_sum=1e4); sc.pp.log1p(B)
sc.pp.calculate_qc_metrics(B, percent_top=None, log1p=False, inplace=True)
B.obs["log_n_genes"] = np.log1p(B.obs["n_genes_by_counts"].astype(float))
sc.pp.regress_out(B, ["log_n_genes"])            # 回归掉深度差异
sc.pp.scale(B, max_value=10)
sc.tl.pca(B, n_comps=50, svd_solver="arpack")
sc.external.pp.harmony_integrate(B, "dataset", basis="X_pca", adjusted_basis="X_pca_harmony", max_iter_harmony=20)
sc.pp.neighbors(B, use_rep="X_pca_harmony", n_neighbors=15, n_pcs=30)
print("=== 加强版（3000 HVG + 回归深度 + Harmony 20 轮）跨数据集连边 %.1f%% ===" % (mix_edges(B)*100))
sc.tl.leiden(B, resolution=0.5, key_added="cl2", flavor="igraph", n_iterations=2, directed=False)
t = pd.crosstab(B.obs.cl2, B.obs.dataset); t["我的占比%"] = (t["Mine"]/t.sum(axis=1)*100).round(1)
print("\n=== 加强版各群组成 ==="); print(t.to_string())
# 未整合的对照（直接 PCA）
sc.pp.neighbors(A.copy() if False else B, use_rep="X_pca", n_neighbors=15, n_pcs=30, key_added="raw_nb")
C2 = B.copy(); sc.pp.neighbors(C2, use_rep="X_pca", n_neighbors=15, n_pcs=30)
print("=== 未整合对照（仅 PCA）跨数据集连边 %.1f%% ===" % (mix_edges(C2)*100))
sc.tl.umap(B)
B.write(out + r"\joint_mine_elife_v2.h5ad")
# 出图
sc.settings.figdir = out + r"\figures_joint"
pal = {"Mine":"#C44E52", "eLife":"#4C72B0"}
fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
um = B.obsm["X_umap"]; ds = B.obs["dataset"].values
for k, ax in enumerate(axes):
    if k < 2:
        for d in ("Mine","eLife"):
            m = ds == d
            ax.scatter(um[m,0], um[m,1], s=0.7, c=pal[d], linewidths=0, rasterized=True, label="%s (n=%d)" % (d, m.sum()))
        ax.legend(markerscale=8, fontsize=9); ax.set_title("联合 UMAP（Harmony 整合后）" if k==0 else "同上（按数据集着色）", fontsize=11)
    else:
        mm = ds == "Mine"
        ax.scatter(um[~mm,0], um[~mm,1], s=0.7, c="#DDDDDD", linewidths=0, rasterized=True, label="eLife")
        ax.scatter(um[mm,0], um[mm,1], s=0.7, c=pal["Mine"], linewidths=0, rasterized=True, label="Mine")
        ax.legend(markerscale=8, fontsize=9); ax.set_title("我的数据（红色）vs eLife（灰）", fontsize=11)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
plt.tight_layout(); plt.savefig(out + r"\figures_joint\joint_umap.png", dpi=150); plt.close()
print("\nfigure saved")
