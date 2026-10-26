import sys, io, os
import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, scanpy as sc, anndata as ad, scipy.sparse as sp
from scipy.stats import rankdata
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 250); sc.settings.verbosity = 0
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
S_G = ['PCNA','MCM6','CDC6','ORC1','ORC4','ORC5','POLA2','PRIM1','FEN1','TYMS','DHFR','HIRA','CHAF1A']
G2M_G = ['CCNB3','AURKA','TOP2A','TPX2','KIF20A','KIF23','RACGAP1','NDC80','INCENP','PRC1','ASPM','Pka-C1','Patronin']
A = sc.read_h5ad(out + r"\adata_soupx_hvg.h5ad")      # Harmony 坐标 + HVG
Araw = sc.read_h5ad(out + r"\adata_soupx.h5ad")
Araw.X = Araw.layers["counts"].copy(); sc.pp.normalize_total(Araw, target_sum=1e4); sc.pp.log1p(Araw)
sym2i = {}
for i, s in zip(Araw.var_names, Araw.var["sym"]): sym2i.setdefault(str(s), i)
sc.tl.score_genes_cell_cycle(Araw, s_genes=[sym2i[g] for g in S_G], g2m_genes=[sym2i[g] for g in G2M_G])
H = Araw[Araw.obs.vendor_cluster.astype(str).isin(["0","1","2"])].copy()
print("血细胞 subset:", H.shape)
print("\n=== 各原始群细胞周期（校正后）===")
print(pd.crosstab(H.obs.vendor_cluster.astype(str), H.obs["phase"], normalize="index").mul(100).round(1).to_string())
print(H.obs.groupby(H.obs.vendor_cluster.astype(str))[["S_score","G2M_score"]].mean().round(3).to_string())
# 用校正后的 harmony 坐标做子集重算
sub_idx = [list(A.obs_names).index(x) for x in H.obs_names]
H.obsm["X_pca_harmony"] = A.obsm["X_pca_harmony"][sub_idx]
sc.pp.neighbors(H, use_rep="X_pca_harmony", n_neighbors=15, n_pcs=30)
sc.tl.paga(H, groups="vendor_cluster")
ct = pd.DataFrame(H.uns["paga"]["connectivities"].todense(), index=H.obs.vendor_cluster.cat.categories, columns=H.obs.vendor_cluster.cat.categories)
print("\n=== 血细胞 PAGA（校正后）==="); print(ct.round(3).to_string())
sc.tl.diffmap(H, n_comps=15)
sc_ = H.obs.groupby(H.obs.vendor_cluster.astype(str))[["S_score","G2M_score"]].mean(); sc_["tot"]=sc_.sum(axis=1)
rootc = sc_["tot"].idxmax(); idx = np.where(H.obs.vendor_cluster.astype(str).values==rootc)[0]
H.uns["iroot"] = int(idx[np.argmax(H.obs["G2M_score"].values[idx])]); sc.tl.dpt(H)
print("\n根群 = C%s | 各群 pseudotime 中位:" % rootc)
print(H.obs.groupby(H.obs.vendor_cluster.astype(str))["dpt_pseudotime"].median().round(3).to_string())
for g in ["blank","treated"]:
    v = H.obs.loc[H.obs.group==g,"dpt_pseudotime"]
    print("  %-8s 中位=%.3f | 晚期(>0.5)=%.1f%%" % (g, v.median(), (v>0.5).mean()*100))
X = H.X; pt = H.obs["dpt_pseudotime"].values
rng = np.random.default_rng(0); si = rng.choice(len(pt), size=min(4000,len(pt)), replace=False)
Xs, pts = X[si], pt[si]
ry = rankdata(pts); ry = (ry-ry.mean())/(ry.std()+1e-12)
mean_ = np.asarray(Xs.mean(0)).ravel(); sq_ = np.asarray(Xs.multiply(Xs).mean(0)).ravel(); var = sq_ - mean_**2
hv = np.argsort(-var)[:4000]
Xs = Xs[:, hv].toarray().astype(np.float32)
rx = rankdata(Xs, axis=0); rx = (rx-rx.mean(0))/(rx.std(0)+1e-12)
rho = (rx*ry[:,None]).mean(0)
res = pd.DataFrame({"gene": H.var["sym"].values[hv].astype(str), "rho": rho}).dropna().sort_values("rho")
res.to_csv(out + r"\genes_vs_pseudotime_soupx.csv", index=False)
print("\n=== 校正后：随拟时序下降（祖细胞特征）top12 ===")
print(", ".join("%s(%.2f)" % (r.gene, r.rho) for _, r in res.head(12).iterrows()))
print("=== 校正后：随拟时序上升（成熟特征）top12 ===")
print(", ".join("%s(%.2f)" % (r.gene, r.rho) for _, r in res.tail(12).iloc[::-1].iterrows()))
H.write(out + r"\adata_hemocyte_soupx.h5ad")
print("\nsaved adata_hemocyte_soupx.h5ad")

