# -*- coding: utf-8 -*-
"""155_C0_markers.py - C0 的 top marker（含符号、检出率、log2FC），并做功能归类"""
import io, sys, os, re
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"
A = sc.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))
A.X = A.layers["counts"].copy(); sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
ct = A.obs["vendor_cluster"].astype(str)
sub = A[ct.isin(["0","1","2","3"]).values].copy()          # 只在血细胞四群里比较
sub.obs["ct"] = ct[ct.isin(["0","1","2","3"])].values
sym = sub.var["sym"].astype(str).values
sc.tl.rank_genes_groups(sub, "ct", groups=["0"], reference="rest", method="wilcoxon", pts=True)
d = sc.get.rank_genes_groups_df(sub, group="0")
d["symbol"] = d["names"].map(dict(zip(sub.var_names.astype(str), sym)))
d["pct_diff"] = (d.pct_nz_group - d.pct_nz_reference) * 100
d = d.sort_values("scores", ascending=False)
print("=== C0 在血细胞四群中的 top 40 marker ===")
print(d.head(40)[["symbol","names","logfoldchanges","pct_nz_group","pct_nz_reference","pct_diff","pvals_adj"]].round(3).to_string(index=False))
d.to_csv(os.path.join(OUT, "C0_top_markers.csv"), index=False)

# 功能粗分类（按基因名关键词）
def cls(g):
    g = str(g)
    if re.match(r"^(Rp[LS]|rpl-|rps-|RpLP|eef|eIF|RACK)", g, re.I): return "翻译/核糖体"
    if re.match(r"^(Act|Tub|Myo|Mhc|Cap|Cof|Arp|Fil|Parv|Zasp|Lasp|Talin|Vinc)", g, re.I): return "细胞骨架/黏附"
    if re.match(r"^(PEN|crustin|ALF|A2M|PPAF|TGM|SSPO|hemocytin|Lectin|CTL|Gal)", g, re.I): return "免疫/抗菌/凝血"
    if re.match(r"^(Cts|Prss|Try|Cpa|Chia|Amy|LCP|CTSC|legumain|LGMN)", g, re.I): return "蛋白酶/消化"
    if re.match(r"^(Nlrp|Casp|TNF|TLR|MyD|Rel|IKK|JNK|MAPK|STAT|Dorsal)", g, re.I): return "信号/免疫通路"
    if re.match(r"^(SLC|ABC|ATP|COX|ND|CYTB|SDH|MDH|GAPDH|PGK|PK|HK|LDH)", g, re.I): return "代谢/转运"
    if re.match(r"^(Hsp|hsc|DNAJ|GST|SOD|CAT|GPX|PRDX|Trx)", g, re.I): return "应激/抗氧化"
    if re.match(r"^(Hist|H1|H2|H3|H4|Dsp|HMGB)", g, re.I): return "染色质/核"
    return "其他/未分类"
top = d.head(200).copy(); top["class"] = top["symbol"].map(cls)
print("\n=== top 200 的功能分布 ===")
print(top["class"].value_counts().to_string())
print("\n=== 每类前 6 个基因 ===")
for k, g in top.groupby("class"):
    print("  %-12s %s" % (k, ", ".join(g["symbol"].astype(str).head(6))))
