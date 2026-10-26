# -*- coding: utf-8 -*-
"""154_elife_scoregenes_percluster.py - score_genes（背景校正）在 10 个群上的均值，检验是否 GH1 普遍偏高"""
import io, sys, os
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"
EMAP = {"elife-80127-supp3-v2.csv":"PH1","elife-80127-supp4-v2.csv":"GH2",
        "elife-80127-supp5-v2.csv":"MH2","elife-80127-supp6-v2.csv":"PH2",
        "elife-80127-supp7-v2.csv":"GH1","elife-80127-supp8-v2.csv":"MH1"}
OUR = {"0":"MH1","1":"GH","2":"PH","3":"MH2","4":"Epi","5":"F","6":"R","7":"B","8":"Neu8","9":"Neu9"}
A = sc.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))
A.X = A.layers["counts"].copy(); sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
ids = list(A.var_names.astype(str)); id2i = {g: i for i, g in enumerate(ids)}
sets = {}
for f, name in EMAP.items():
    d = pd.read_csv(os.path.join(OUT, "elife_ref", f))
    if "Log2FC" in d.columns: d = d.sort_values("Log2FC", ascending=False)
    sets[name] = [g for g in d["Gene ID"].astype(str) if g in id2i][:100]
for k, gl in sets.items():
    sc.tl.score_genes(A, gl, ctrl_size=50, score_name="s_" + k, random_state=0)
cols = ["s_" + k for k in sets]
cl = A.obs["vendor_cluster"].astype(str).values
tab = pd.DataFrame({k: [A.obs.loc[cl == c, "s_" + k].mean() for c in OUR] for k in sets}, index=list(OUR.values()))
tab["n"] = [(cl == c).sum() for c in OUR]
tab["best"] = tab[list(sets)].idxmax(axis=1)
print("=== score_genes（背景校正，ctrl=50）在各群上的均值 ===")
print(tab.round(3).to_string())
print("\n各集合的群间极差（最大值-最小值）:")
for k in sets:
    print("  %-4s range %.3f  (max %s %.3f, min %s %.3f)" % (k, tab[k].max()-tab[k].min(),
          tab[k].idxmax(), tab[k].max(), tab[k].idxmin(), tab[k].min()))
tab.to_csv(os.path.join(OUT, "elife_scoregenes_by_cluster.csv"))
