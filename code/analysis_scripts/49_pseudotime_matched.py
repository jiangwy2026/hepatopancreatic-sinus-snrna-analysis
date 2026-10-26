import sys, io, os
import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, scanpy as sc
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 250); sc.settings.verbosity = 0
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
H = sc.read_h5ad(out + r"\adata_hemocyte_traj.h5ad")
H.X = H.layers["counts"].copy(); sc.pp.normalize_total(H, target_sum=1e4); sc.pp.log1p(H)
sym2i = {}
for i, s in zip(H.var_names, H.var["sym"].astype(str)): sym2i.setdefault(s, i)
c2 = H[H.obs.vendor_cluster.astype(str) == "2"].copy()
# 1) C2 内部的 C0 样信号（substructure 检查）
C0MK = ["THBS4","LOC113830563","LOC113830564","LOC113830566","LOC113819543","LOC113805372","LOC113823493"]
ids = [sym2i[g] for g in C0MK if g in sym2i]
c2.obs["C0_score"] = np.asarray(c2[:, ids].X.mean(1)).ravel()
C2MK = ["LOC113823945","ccnd2","Igfbp7","Grik2","slo"]
ids2 = [sym2i[g] for g in C2MK if g in sym2i]
c2.obs["C2_score"] = np.asarray(c2[:, ids2].X.mean(1)).ravel()
print("=== C2 群内部：组间 substructure 检查 ===")
print(c2.obs.groupby("group")[["C0_score","C2_score","G2M_score","S_score"]].mean().round(3).to_string())
print("\n  检出率>0 的 C0 marker 个数/细胞:")
c2.obs["n_C0mk"] = np.asarray((c2[:, ids].X > 0).sum(1)).ravel()
print(c2.obs.groupby("group")["n_C0mk"].mean().round(2).to_string())
# 2) 拟时序匹配比较：按 pseudotime 分箱后，箱内比较组间
print("\n=== 拟时序匹配后：C2 内部组间差异（按 pseudotime 三分位）===")
c2.obs["pt_bin"] = pd.qcut(c2.obs["dpt_pseudotime"], 3, labels=["早(P1)","中(P2)","晚(P3)"])
for b in ["早(P1)","中(P2)","晚(P3)"]:
    sub = c2[c2.obs.pt_bin == b]
    n_b = (sub.obs.group=="blank").sum(); n_t = (sub.obs.group=="treated").sum()
    if n_b < 20 or n_t < 20: print("  %s: 细胞太少 (blank %d / treated %d)" % (b, n_b, n_t)); continue
    sc.tl.rank_genes_groups(sub, "group", groups=["treated"], reference="blank", method="wilcoxon", pts=True)
    dd = sc.get.rank_genes_groups_df(sub, group="treated")
    dd["sym"] = dd["names"].map(dict(zip(sub.var_names, sub.var["sym"].astype(str))))
    dd["gene"] = dd["sym"].fillna(dd["names"]).astype(str)
    s = dd[(dd.pvals_adj<0.05)&(dd.logfoldchanges>0.5)]
    print("\n  [%s] blank %d / treated %d | 处理组上调 %d 个, top8: %s" % (b, n_b, n_t, len(s), ", ".join(s.sort_values("scores",ascending=False).head(8).gene.tolist())))
    dn = dd[(dd.pvals_adj<0.05)&(dd.logfoldchanges<-0.5)]
    print("      处理组下调 %d 个, top8: %s" % (len(dn), ", ".join(dn.sort_values("scores").head(8).gene.tolist())))
# 3) 只在“早”箱里，且组间一致上调的基因（更可能驱动扩增）
sub = c2[c2.obs.pt_bin == "早(P1)"]
sc.tl.rank_genes_groups(sub, "group", groups=["treated"], reference="blank", method="wilcoxon", pts=True)
d0 = sc.get.rank_genes_groups_df(sub, group="treated")
d0["sym"] = d0["names"].map(dict(zip(sub.var_names, sub.var["sym"].astype(str))))
d0["gene"] = d0["sym"].fillna(d0["names"]).astype(str)
tf = {}
with open(root + r"\all.genes.marker.stat.xls.xls", encoding="utf-8", errors="replace") as fh:
    hd = fh.readline().rstrip("\n").split("\t"); i_s, i_t = hd.index("Symbol"), hd.index("TF_family")
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) > i_t and p[i_t] not in ("-",""): tf[p[i_s]] = p[i_t]
d0["TF"] = d0["gene"].map(tf)
s = d0[(d0.pvals_adj<0.05)&(d0.logfoldchanges>0.5)].sort_values("scores", ascending=False)
print("\n=== 早期祖细胞里处理组上调的基因（候选扩增驱动因子）===")
for _, r in s.head(20).iterrows():
    print("   %-22s log2FC=%+.2f q=%.1e %s" % (r.gene, r.logfoldchanges, r.pvals_adj, ("[TF:%s]"%r.TF) if pd.notna(r.TF) else ""))
s.to_csv(out + r"\C2_early_treated_up.csv", index=False)
c2.write(out + r"\adata_C2_within.h5ad")
