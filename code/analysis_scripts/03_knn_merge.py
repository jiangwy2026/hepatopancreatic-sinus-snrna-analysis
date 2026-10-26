import sys, io, os
import numpy as np, pandas as pd, scanpy as sc
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
fig = out + r"\figures"; os.makedirs(fig, exist_ok=True)
B = sc.read_h5ad(out + r"\adata_hvg_processed.h5ad")
lab = B.obs["vendor_cluster"].astype(str).values
order = sorted(set(lab), key=int); name = ["C" + o for o in order]
codes = np.array([order.index(x) for x in lab]); k = len(order)
C = B.obsp["connectivities"].tocoo()
R = np.zeros((k, k)); np.add.at(R, (codes[C.row], codes[C.col]), C.data); R = (R + R.T) / 2
rown = R / R.sum(1, keepdims=True)
D = lambda M, nd=1: pd.DataFrame(np.round(M, nd), index=name, columns=name)
print("=== 每群 kNN 连通性去向（行归一化 %，对角线=自身）===")
print(D(rown * 100, 1).to_string())
pairs = sorted(((order[i], order[j], (rown[i,j]+rown[j,i])/2*100) for i in range(k) for j in range(i+1,k)), key=lambda x: -x[2])
print("\n=== 最该合并的群对（双向平均连通占比）===")
for a,b,v in pairs[:12]: print("   C%-2s - C%-2s : %.2f%%" % (a,b,v))
sizes = np.bincount(codes, minlength=k)
exp = np.outer(sizes, sizes) / (sizes.sum()**2); obs = R / R.sum()
enr = np.where(exp > 0, obs/exp, np.nan)
print("\n=== 群间连接富集倍数 obs/exp（>3 说明高度互连 = 合并候选）===")
print(D(enr, 2).to_string())
D(enr, 3).to_csv(out + r"\cluster_knn_enrichment.csv")
ct = pd.crosstab(B.obs["vendor_cluster"].astype(str), B.obs["leiden_0.3"].astype(str))
print("\n=== vendor cluster x leiden(0.3) ===")
print(ct.to_string()); ct.to_csv(out + r"\crosstab_vendor_vs_leiden03.csv")
sc.settings.figdir = fig
sc.pl.umap(B, color="vendor_cluster", legend_loc="on data", frameon=True, title="Vendor clusters (0-9)", save="_vendor_cluster.png", show=False, palette="tab10")
sc.pl.umap(B, color="group", frameon=True, title="Group", save="_group.png", show=False)
sc.pl.umap(B, color="sample", frameon=True, title="Sample", save="_sample.png", show=False)
print("\nfigures:", os.listdir(fig))
