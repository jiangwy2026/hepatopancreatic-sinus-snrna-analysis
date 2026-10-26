import sys, io, os
import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, scanpy as sc
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 250)
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"; fig = out + r"\figures_sub"
os.makedirs(fig, exist_ok=True)
Hh = sc.read_h5ad(out + r"\adata_hemocyte_harmony.h5ad")
sub = Hh.obs["sub_0.3"].astype(str)
umap = pd.DataFrame(Hh.obsm["X_umap"] if "X_umap" in Hh.obsm else np.zeros((Hh.n_obs,2)), index=Hh.obs_names)

# 在“未 scaled、log-normalized 全基因”对象上算 marker
A = sc.read_h5ad(out + r"\adata_raw.h5ad")
A = A[A.obs_names.isin(set(Hh.obs_names))].copy()
A.X = A.layers["counts"].copy()
sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
A.obs["sub"] = sub.loc[A.obs_names].values
sc.tl.rank_genes_groups(A, "sub", method="wilcoxon", pts=True)
mk = pd.concat([sc.get.rank_genes_groups_df(A, group=g).assign(sub=g) for g in sorted(A.obs["sub"].unique())], ignore_index=True)
mk = mk.merge(A.var.reset_index()[["index","symbol"]].rename(columns={"index":"names"}), on="names", how="left")
mk["gene"] = mk["symbol"].fillna(mk["names"])
mk.to_csv(out + r"\subcluster_markers_final.csv", index=False)
sig = mk[(mk.pvals_adj < 0.05) & (mk.logfoldchanges > 0.5)].copy()
print("=== 每亚群显著 marker 数 ==="); print(sig.groupby("sub").size().to_string())
for g in sorted(A.obs["sub"].unique()):
    d = sig[sig["sub"] == g].sort_values("scores", ascending=False).head(15)
    print("\n【亚群 %s】 n=%d  (原分群: %s)" % (g, (A.obs["sub"]==g).sum(),
          dict(pd.crosstab(A.obs["sub"], Hh.obs["vendor_cluster"]).loc[g])))
    for _, r in d.iterrows():
        print("   %-22s log2FC=%+.2f  padj=%.1e  %2.0f%%/%2.0f%%" % (r.gene, r.logfoldchanges, r.pvals_adj, r.pct_nz_group*100, r.pct_nz_reference*100))
sig.sort_values(["sub","scores"],ascending=[True,False]).groupby("sub").head(40).to_csv(out+r"\subcluster_top40_final.csv", index=False, encoding="utf-8-sig")

# UMAP 图（用 harmony 邻接图重算）
sc.pp.neighbors(Hh, use_rep="X_pca_harmony", n_neighbors=15, n_pcs=30, key_added="hn2")
sc.tl.umap(Hh, neighbors_key="hn2")
sc.settings.figdir = fig
sc.pl.umap(Hh, color="sub_0.3", legend_loc="on data", title="Hemocyte subtypes (Harmony, res=0.3)", save="_sub_cluster.png", show=False, palette="tab10")
sc.pl.umap(Hh, color="vendor_cluster", title="Original clusters (0/1/2)", save="_vendor.png", show=False)
sc.pl.umap(Hh, color="sample", title="Sample (Harmony-integrated)", save="_sample.png", show=False)
Hh.write(out + r"\adata_hemocyte_harmony.h5ad")
print("\nfigures:", sorted(os.listdir(fig)))
