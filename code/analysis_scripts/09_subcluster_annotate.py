import sys, io
import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, scanpy as sc
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 250)
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
H = sc.read_h5ad(out + r"\adata_hemocyte_harmony.h5ad")
q = pd.read_csv(out + r"\qc_per_cell_raw.csv").set_index("Cells")
H.obs["nCount"] = q.loc[H.obs_names, "nCount"].values
lab = "sub_0.3"
print("=== 亚群 x 样本（整合后，应充分混合）===")
print(pd.crosstab(H.obs[lab], H.obs["sample"]).to_string())
t = pd.crosstab(H.obs[lab], H.obs["group"])
t["pct_blank"]=(t["blank"]/t["blank"].sum()*100).round(1); t["pct_treated"]=(t["treated"]/t["treated"].sum()*100).round(1)
t["fold"]=(t["treated"]/t["blank"]).round(2)
print("\n=== 亚群 x 组 ==="); print(t.to_string())
print("\n=== 亚群 QC 中位数 ==="); print(H.obs.groupby(lab)[["nCount","n_genes_by_counts","pct_counts_mt"]].median().round(1).to_string())
print("\n=== 亚群 x 原 血细胞分群 ==="); print(pd.crosstab(H.obs[lab], H.obs.vendor_cluster).to_string())

sc.tl.rank_genes_groups(H, lab, method="wilcoxon", pts=True, use_raw=False)
mk = pd.concat([sc.get.rank_genes_groups_df(H, group=g).assign(sub=g) for g in H.obs[lab].cat.categories], ignore_index=True)
mk = mk.merge(H.var.reset_index()[["index","symbol"]].rename(columns={"index":"names"}), on="names", how="left")
mk["gene"] = mk["symbol"].fillna(mk["names"])
mk.to_csv(out + r"\subcluster_markers_harmony.csv", index=False)
sig = mk[(mk.pvals_adj<0.05) & (mk.logfoldchanges>0.5)]
print("\n=== 每亚群显著 marker 数 ==="); print(sig.groupby("sub").size().to_string())
print("\n=== 每亚群 top12 marker（含检出率）===")
for g in H.obs[lab].cat.categories:
    d = sig[sig.sub==g].sort_values("scores", ascending=False).head(12)
    print("\n亚群 %s (n=%d):" % (g, (H.obs[lab]==g).sum()))
    for _, r in d.iterrows():
        print("   %-20s log2FC=%+.2f  padj=%.1e  本群%.0f%% / 其余%.0f%%" % (r.gene, r.logfoldchanges, r.pvals_adj, r.pct_nz_group*100, r.pct_nz_reference*100))
sig.sort_values(["sub","scores"],ascending=[True,False]).groupby("sub").head(40).to_csv(out+r"\subcluster_top40_harmony.csv", index=False, encoding="utf-8-sig")
print("\nwrote subcluster_markers_harmony.csv / subcluster_top40_harmony.csv")
