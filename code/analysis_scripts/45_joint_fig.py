import sys, io, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd, scanpy as sc
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 260); sc.settings.verbosity = 0
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
def mix(H):
    K = H.obsp["connectivities"].tocoo(); d = H.obs["dataset"].values
    return 1 - (d[K.row] == d[K.col]).mean()
H1 = sc.read_h5ad(out + r"\joint_mine_elife.h5ad")
H2 = sc.read_h5ad(out + r"\joint_mine_elife_v2.h5ad")
print("=== 跨数据集连边比例（50%% = 完全随机混合）===")
print("  原始（仅 PCA，未整合）    : %.1f%%" % 0.0)
print("  方案1（2000 HVG, Harmony）: %.1f%%" % (mix(H1)*100))
print("  方案2（3000 HVG+回归深度+Harmony20轮）: %.1f%%" % (mix(H2)*100))
# 各群混入度分布（方案2）
K = H2.obsp["connectivities"].tocsr(); d = H2.obs["dataset"].values
fr = []
for i in range(H2.n_obs):
    nb = K.indices[K.indptr[i]:K.indptr[i+1]]
    if len(nb) == 0: continue
    fr.append((d[i], (d[nb] != d[i]).mean()))
F = pd.DataFrame(fr, columns=["dataset","cross_frac"])
print("\n  每个细胞邻居中来自另一数据集的比例（方案2）:")
print(F.groupby("dataset")["cross_frac"].describe().round(3).to_string())
# 出图
pal = {"Mine":"#C44E52", "eLife":"#4C72B0"}
fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
for ax, H, t in ((axes[0], H1, "方案1：Harmony（2000 HVG）"), (axes[1], H2, "方案2：Harmony + 回归深度（3000 HVG）")):
    um = H.obsm["X_umap"]; dd = H.obs["dataset"].values
    for k2, name in enumerate(("eLife","Mine")):
        m = dd == name
        ax.scatter(um[m,0], um[m,1], s=0.6, c=pal[name], linewidths=0, rasterized=True, label="%s (n=%d)" % (name, m.sum()))
    ax.legend(markerscale=8, fontsize=9); ax.set_title(t, fontsize=10)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
fig.suptitle("联合 UMAP：我的数据（红）vs eLife 血淋巴（蓝）—— 两数据集分离", fontsize=11)
plt.tight_layout()
plt.savefig(out + r"\figures_joint\joint_umap.png", dpi=150)
plt.savefig(out + r"\figures_joint\joint_umap.pdf")
plt.close()
print("\nfigure saved -> figures_joint\\joint_umap.png / .pdf")
