import sys, io, os
import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, scanpy as sc
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 250)
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
Hs = sc.read_h5ad(out + r"\adata_hemocyte_sub.h5ad")
lab = "leiden_0.5"
print("=== 亚群 x 原分群 ==="); print(pd.crosstab(Hs.obs[lab], Hs.obs.vendor_cluster).to_string())
print("\n=== 亚群 x 样本组（细胞数 / 占比）===")
t = pd.crosstab(Hs.obs[lab], Hs.obs.group)
t["pct_blank"] = (t["blank"]/t["blank"].sum()*100).round(1); t["pct_treated"] = (t["treated"]/t["treated"].sum()*100).round(1)
t["fold"] = (t["treated"]/t["blank"]).round(2); print(t.to_string())
print("\n=== 亚群 x 样本（原始细胞数），检查是否被深度/批次驱动 ===")
print(pd.crosstab(Hs.obs[lab], Hs.obs["sample"]).to_string())
A0 = sc.read_h5ad(out + r"\adata_raw.h5ad"); q = pd.read_csv(out + r"\qc_per_cell_raw.csv").set_index("Cells")
Hs.obs["nCount"] = q.loc[Hs.obs_names, "nCount"].values
print("\n=== 各亚群 QC 中位数 ===")
print(Hs.obs.groupby(lab)[["nCount","n_genes_by_counts","pct_counts_mt"]].median().round(1).to_string())
print("\n=== Wilcoxon marker（padj<0.05 & log2FC>0.5），每亚群显著数 ===")
Hn = Hs.copy(); Hn.X = np.expm1(Hn.X) if False else Hn.X
sc.tl.rank_genes_groups(Hs, lab, method="wilcoxon", pts=True)
mk = pd.concat([sc.get.rank_genes_groups_df(Hs, group=g).assign(sub=g) for g in Hs.obs[lab].cat.categories], ignore_index=True)
mk = mk.merge(Hs.var.reset_index()[["index","symbol"]].rename(columns={"index":"names"}), on="names", how="left")
mk["gene"] = mk["symbol"].fillna(mk["names"])
sig = mk[(mk.pvals_adj<0.05) & (mk.logfoldchanges>0.5)]
print(sig.groupby("sub").size().to_string())
for g in Hs.obs[lab].cat.categories:
    d = sig[sig.sub==g].sort_values("scores", ascending=False).head(12)
    print("\n  亚群 %-2s: %s" % (g, ", ".join(d.gene.astype(str))))
mk.to_csv(out + r"\subcluster_markers_all.csv", index=False)
sig.sort_values(["sub","scores"], ascending=[True,False]).groupby("sub").head(30).to_csv(out + r"\subcluster_top30_markers.csv", index=False, encoding="utf-8-sig")
print("\nwrote subcluster_markers_all.csv / subcluster_top30_markers.csv")
