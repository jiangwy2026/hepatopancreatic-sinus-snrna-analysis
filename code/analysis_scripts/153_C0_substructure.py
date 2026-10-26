# -*- coding: utf-8 -*-
"""153_C0_substructure.py - C0 内部是否有亚结构？各亚结构分别更像哪个 eLife 亚型？"""
import io, sys, os
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"
EMAP = {"elife-80127-supp3-v2.csv":"PH1","elife-80127-supp4-v2.csv":"GH2",
        "elife-80127-supp5-v2.csv":"MH2","elife-80127-supp6-v2.csv":"PH2",
        "elife-80127-supp7-v2.csv":"GH1","elife-80127-supp8-v2.csv":"MH1"}
A = sc.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))
A.X = A.layers["counts"].copy(); sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
ids = list(A.var_names.astype(str)); id2i = {g: i for i, g in enumerate(ids)}
C = A[A.obs["vendor_cluster"].astype(str) == "0"].copy()
print("C0 细胞数:", C.n_obs, "| blank/treated:", C.obs["sample"].value_counts().to_dict())
sc.pp.highly_variable_genes(C, n_top_genes=1500); Ch = C[:, C.var.highly_variable].copy()
sc.pp.scale(Ch, max_value=10); sc.tl.pca(Ch, n_comps=30, svd_solver="arpack")
sc.pp.neighbors(Ch, n_neighbors=15, use_rep="X_pca")
for res in (0.2, 0.4):
    sc.tl.leiden(Ch, resolution=res, flavor="igraph", n_iterations=2, directed=False, key_added="c0_%.1f" % res)
    print("\nres %.1f -> %d 个亚群" % (res, Ch.obs["c0_%.1f" % res].nunique()))
key = "c0_0.4"
print(Ch.obs.groupby([key, "sample"], observed=True).size().unstack(fill_value=0).to_string())
sets = {}
for f, name in EMAP.items():
    d = pd.read_csv(os.path.join(OUT, "elife_ref", f))
    if "Log2FC" in d.columns: d = d.sort_values("Log2FC", ascending=False)
    sets[name] = [g for g in d["Gene ID"].astype(str) if g in id2i][:150]
S = pd.DataFrame(index=C.obs_names)
for k, gl in sets.items():
    sc.tl.score_genes(C, gl, ctrl_size=100, score_name="s_" + k, random_state=0)
    S[k] = C.obs["s_" + k].values
S[key] = Ch.obs[key].values
tab = S.groupby(key).mean(numeric_only=True)
tab["n"] = S.groupby(key).size()
tab["best"] = tab[list(sets)].idxmax(axis=1)
print("\n=== C0 各亚群的 eLife 打分均值 ===")
print(tab.round(3).to_string())
tab.to_csv(os.path.join(OUT, "C0_substructure_elife.csv"))
