# -*- coding: utf-8 -*-
"""156_C0_ribosome_check.py - C0 是不是"高核糖体/高 RNA"状态？并查 C0 特异 LOC 基因的注释"""
import io, sys, os, re, gzip
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"
A = sc.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))
X = A.layers["counts"].tocsr()
cl = A.obs["vendor_cluster"].astype(str).values
sym = A.var["sym"].astype(str).values
rib = np.array([bool(re.match(r"^(Rp[LS]|rpl-|rps-|RpLP|RPS|RPL)", s, re.I)) for s in sym])
print("核糖体基因数:", rib.sum())
tot = np.asarray(X.sum(1)).ravel()
ng = np.asarray((X > 0).sum(1)).ravel()
ribc = np.asarray(X[:, rib].sum(1)).ravel()
print("\n=== 各群：UMI、基因数、核糖体占比 ===")
rows = []
for c in ["0","1","2","3","4","5","6","7","8","9"]:
    m = cl == c
    rows.append({"cluster": c, "n": m.sum(), "UMI_med": np.median(tot[m]), "n_genes_med": np.median(ng[m]),
                 "ribo_pct_med": round(100*np.median(ribc[m]/np.maximum(tot[m],1)), 1),
                 "ribo_umi_med": np.median(ribc[m])})
print(pd.DataFrame(rows).to_string(index=False))

# 用逻辑回归/相关看"核糖体占比"与群的关系：是否只由 UMI 决定
print("\n=== 按 UMI 分层后，C0 与其他群的核糖体占比 ===")
q = pd.qcut(tot, 5, labels=["Q1","Q2","Q3","Q4","Q5"])
df = pd.DataFrame({"c": cl, "q": q, "ribo": 100*ribc/np.maximum(tot,1), "umi": tot})
print(df[df.c.isin(["0","1","2","3"])].groupby(["q","c"], observed=True)["ribo"].median().unstack().round(1).to_string())

# C0 特异 LOC 基因的注释
locs = ["LOC113819543","LOC113823493","LOC113830563","LOC113830564","LOC113830566","LOC113815249","LOC113819837","LOC113827561","LOC113805372","LOC113830566"]
gtf = r"E:\workspace\genome2\genome\GENOME_GTF\GCF_003789085.1.gtf.gz"
print("\n=== GTF 中这些 LOC 基因的注释 ===")
found = {}
with gzip.open(gtf, "rt") as f:
    for line in f:
        if line.startswith("#"): continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] != "gene": continue
        m = re.search(r'gene_id "([^"]+)"', p[8])
        if m and m.group(1) in locs:
            found[m.group(1)] = p[8].strip()
for g in locs:
    print("  %-14s %s" % (g, found.get(g, "(GTF 未找到)")))
# eLife 补充表里的 Description
print("\n=== eLife 补充表中的 Description（若有）===")
import glob
desc = {}
for f in glob.glob(os.path.join(OUT, "elife_ref", "*.csv")):
    d = pd.read_csv(f)
    if "Gene ID" in d.columns and "Description" in d.columns:
        for r in d.itertuples():
            if str(r._2) in locs:
                desc[str(r._2)] = str(r.Description)
for g in locs:
    if g in desc: print("  %-14s %s" % (g, desc[g][:120]))
