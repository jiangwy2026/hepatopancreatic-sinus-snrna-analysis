import sys, io, os
import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 250); sc.settings.verbosity = 0
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
A = sc.read_h5ad(out + r"\adata_soupx_hvg.h5ad")
S = sc.read_h5ad(out + r"\adata_soupx.h5ad")          # 校正后 counts（全基因）
S.X = S.layers["counts"].copy(); sc.pp.normalize_total(S, target_sum=1e4); sc.pp.log1p(S)
H = S[S.obs.vendor_cluster.astype(str).isin(["0","1","2"])].copy()
print("血细胞 subset:", H.shape)
idx = [list(A.obs_names).index(x) for x in H.obs_names]
H.obsm["X_pca_harmony"] = A.obsm["X_pca_harmony"][idx]
sc.pp.neighbors(H, use_rep="X_pca_harmony", n_neighbors=15, n_pcs=30)
print("\n=== 校正后：血细胞亚型多分辨率（Harmony）===")
for r in [0.1,0.2,0.3,0.4,0.5,0.6,0.8]:
    k = "h_%.1f" % r
    sc.tl.leiden(H, resolution=r, key_added=k, flavor="igraph", n_iterations=2, directed=False)
    s = H.obs[k].value_counts()
    print("  res %.1f -> %2d 群 (min %3d, max %4d)" % (r, len(s), s.min(), s.max()))
lab = "h_0.3"
print("\n=== 亚型 x 原始分群 ==="); print(pd.crosstab(H.obs[lab], H.obs.vendor_cluster.astype(str)).to_string())
print("\n=== 亚型 x 样本（是否被批次驱动）===")
t = pd.crosstab(H.obs[lab], H.obs["group"]); t["pct_blank"]=(t["blank"]/t["blank"].sum()*100).round(1); t["pct_treated"]=(t["treated"]/t["treated"].sum()*100).round(1); t["fold"]=(t["treated"]/t["blank"]).round(2)
print(t.to_string())
sc.tl.rank_genes_groups(H, lab, method="wilcoxon", pts=True)
mk = pd.concat([sc.get.rank_genes_groups_df(H, group=g).assign(sub=g) for g in H.obs[lab].cat.categories], ignore_index=True)
mk = mk.merge(H.var.reset_index()[["index","sym"]].rename(columns={"index":"names"}), on="names", how="left")
mk["gene"] = mk["sym"].fillna(mk["names"]).astype(str)
sig = mk[(mk.pvals_adj<0.05) & (mk.logfoldchanges>0.5)]
print("\n=== 每亚型显著 marker 数 ==="); print(sig.groupby("sub").size().to_string())
for g in H.obs[lab].cat.categories:
    d = sig[sig["sub"]==g].sort_values("scores", ascending=False).head(12)
    print("\n亚型 %-2s (n=%4d): %s" % (g, (H.obs[lab]==g).sum(), ", ".join(d.gene.tolist())))
mk.to_csv(out + r"\hemocyte_soupx_subtype_markers.csv", index=False)
H.write(out + r"\adata_hemocyte_soupx_sub.h5ad")
print("\nsaved")
