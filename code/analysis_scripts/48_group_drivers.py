import sys, io, os
import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, scanpy as sc
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 250); sc.settings.verbosity = 0
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
H = sc.read_h5ad(out + r"\adata_hemocyte_traj.h5ad")
H.X = H.layers["counts"].copy(); sc.pp.normalize_total(H, target_sum=1e4); sc.pp.log1p(H)
print("=== 一、增殖活性：组间比较 ===")
print("\n[全体血细胞] 细胞周期相位分布（%）")
print(pd.crosstab(H.obs["group"], H.obs["phase"], normalize="index").mul(100).round(1).to_string())
print("\n[仅 C2 增殖群] 细胞周期相位分布（%）")
c2 = H[H.obs.vendor_cluster.astype(str) == "2"].copy()
print(pd.crosstab(c2.obs["group"], c2.obs["phase"], normalize="index").mul(100).round(1).to_string())
print("\n[C2] S/G2M 打分均值")
print(c2.obs.groupby("group")[["S_score","G2M_score"]].mean().round(4).to_string())
print("\n[C2] 拟时序分布")
print(c2.obs.groupby("group")["dpt_pseudotime"].describe().round(3).to_string())
print("\n[各群] 拟时序中位数 按组")
print(H.obs.pivot_table(index=H.obs.vendor_cluster.astype(str), columns="group", values="dpt_pseudotime", aggfunc="median").round(3).to_string())
print("\n=== 二、C2 内部：对照 vs 处理 差异基因 ===")
c2.obs["grp2"] = c2.obs["group"].astype(str)
sc.tl.rank_genes_groups(c2, "grp2", groups=["treated"], reference="blank", method="wilcoxon", pts=True)
d = sc.get.rank_genes_groups_df(c2, group="treated")
d["sym"] = d["names"].map(dict(zip(c2.var_names, c2.var["sym"].astype(str))))
d["gene"] = d["sym"].fillna(d["names"]).astype(str)
tf = {}
with open(root + r"\all.genes.marker.stat.xls.xls", encoding="utf-8", errors="replace") as fh:
    hd = fh.readline().rstrip("\n").split("\t"); i_s, i_t = hd.index("Symbol"), hd.index("TF_family")
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) > i_t and p[i_t] not in ("-", ""): tf[p[i_s]] = p[i_t]
d["TF"] = d["gene"].map(tf)
sig = d[(d.pvals_adj < 0.05) & (d.logfoldchanges.abs() > 0.5)].copy()
print("  显著差异基因数: %d" % len(sig))
print("\n  --- 处理组【上调】top15 ---")
for _, r in sig.sort_values("scores", ascending=False).head(15).iterrows():
    print("   %-22s log2FC=%+.2f q=%.1e %s" % (r.gene, r.logfoldchanges, r.pvals_adj, ("[TF:%s]" % r.TF) if pd.notna(r.TF) else ""))
print("\n  --- 处理组【下调】top10 ---")
for _, r in sig.sort_values("scores").head(10).iterrows():
    print("   %-22s log2FC=%+.2f q=%.1e %s" % (r.gene, r.logfoldchanges, r.pvals_adj, ("[TF:%s]" % r.TF) if pd.notna(r.TF) else ""))
sig.to_csv(out + r"\C2_treated_vs_blank_DE.csv", index=False)
d.to_csv(out + r"\C2_treated_vs_blank_all.csv", index=False)
print("\n=== 三、全体血细胞：对照 vs 处理（成分无关的基因变化）===")
sc.tl.rank_genes_groups(H, "group", groups=["treated"], reference="blank", method="wilcoxon", pts=True)
d2 = sc.get.rank_genes_groups_df(H, group="treated")
d2["sym"] = d2["names"].map(dict(zip(H.var_names, H.var["sym"].astype(str))))
d2["gene"] = d2["sym"].fillna(d2["names"]).astype(str); d2["TF"] = d2["gene"].map(tf)
s2 = d2[(d2.pvals_adj<0.05) & (d2.logfoldchanges.abs()>0.5)]
print("  显著: %d | 上调 top12:" % len(s2))
for _, r in s2.sort_values("scores", ascending=False).head(12).iterrows():
    print("   %-22s log2FC=%+.2f %s" % (r.gene, r.logfoldchanges, ("[TF:%s]" % r.TF) if pd.notna(r.TF) else ""))
s2.to_csv(out + r"\hemocyte_treated_vs_blank_DE.csv", index=False)
print("\n=== 四、增殖相关基因的表达（组间）===")
CC = ["PCNA","MCM6","CDC6","ORC1","POLA2","PRIM1","FEN1","TYMS","CCNB3","AURKA","TOP2A","TPX2","KIF20A","KIF23","RACGAP1","ASPM","Pka-C1","ccnd2","Patronin"]
sym2i = {}
for i, s in zip(H.var_names, H.var["sym"].astype(str)): sym2i.setdefault(s, i)
rows = []
for g in CC:
    if g not in sym2i: continue
    i = sym2i[g]
    for grp in ("blank","treated"):
        m = (H.obs["group"]==grp).values & (H.obs.vendor_cluster.astype(str)=="2").values
        v = np.asarray(H[m][:, i].X.todense()).ravel()
        rows.append({"gene":g,"group":grp,"mean":round(float(v.mean()),3),"det%":round(float((v>0).mean()*100),1)})
T = pd.DataFrame(rows).pivot(index="gene", columns="group", values=["mean","det%"])
print(T.round(2).to_string())
