# -*- coding: utf-8 -*-
"""
123_syncytium_probe.py  (Part 1)
19,795 (厂商 filtered) -> 9,658 (本项目) 之间丢掉了什么？特别是高 UMI 的"疑似多核"细胞还在不在？
"""
import io, sys, os
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"

raw = ad.read_h5ad(os.path.join(OUT, "adata_raw.h5ad"))       # 19,795
sx  = ad.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"), backed="r")  # 9,658
kept = set(map(str, sx.obs_names))
raw.obs["kept"] = [str(x) in kept for x in raw.obs_names]
print("厂商 filtered 细胞 %d | 本项目保留 %d | 丢弃 %d (%.1f%%)" % (
    raw.n_obs, raw.obs["kept"].sum(), (~raw.obs["kept"]).sum(), 100 * (~raw.obs["kept"]).mean()))
print()
for lab, m in (("保留", raw.obs["kept"].values), ("丢弃", ~raw.obs["kept"].values)):
    o = raw.obs[m]
    print("%s 细胞 n=%d | total_counts 中位 %.0f 均值 %.0f p99 %.0f 最大 %.0f | n_genes 中位 %.0f | pct_mt 中位 %.1f" % (
        lab, m.sum(), o.total_counts.median(), o.total_counts.mean(), o.total_counts.quantile(.99),
        o.total_counts.max(), o.n_genes_by_counts.median(), o.pct_counts_mt.median()))
print()
# 高 UMI 细胞是否被丢
for thr in (5000, 10000, 20000, 50000):
    m = raw.obs.total_counts.values > thr
    if m.sum() == 0:
        print("UMI > %6d : 0 个" % thr); continue
    print("UMI > %6d : %4d 个 | 其中保留 %4d 个 (%.0f%%) | 丢弃 %4d 个" % (
        thr, m.sum(), (m & raw.obs["kept"].values).sum(), 100 * (m & raw.obs["kept"].values).sum() / m.sum(),
        (m & ~raw.obs["kept"].values).sum()))
print()
# 丢弃细胞的特征（为什么被丢）
d = raw.obs[~raw.obs["kept"].values]
print("丢弃细胞的原因分布：")
print("  n_genes < 200  : %5d (%.0f%%)" % ((d.n_genes_by_counts < 200).sum(), 100 * (d.n_genes_by_counts < 200).mean()))
print("  pct_mt  > 20%%  : %5d (%.0f%%)" % ((d.pct_counts_mt > 20).sum(), 100 * (d.pct_counts_mt > 20).mean()))
print("  total < 500    : %5d (%.0f%%)" % ((d.total_counts < 500).sum(), 100 * (d.total_counts < 500).mean()))
print("  以上都不满足   : %5d (%.0f%%)  <- 既不是低质量也不是高线粒体，需要查明" % (
    ((d.n_genes_by_counts >= 200) & (d.pct_counts_mt <= 20) & (d.total_counts >= 500)).sum(),
    100 * ((d.n_genes_by_counts >= 200) & (d.pct_counts_mt <= 20) & (d.total_counts >= 500)).mean()))
print()
# 保留细胞里的最高 UMI 长什么样
k = raw.obs[raw.obs["kept"].values].sort_values("total_counts", ascending=False).head(12)
print("本项目保留细胞中 UMI 最高的 12 个：")
print(k[["sample", "total_counts", "n_genes_by_counts", "pct_counts_mt"]].to_string())
