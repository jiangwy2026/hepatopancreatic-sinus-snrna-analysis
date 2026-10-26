# -*- coding: utf-8 -*-
"""
126_dropped_clusters.py - 被 QC 删除的 10,137 个细胞是什么？里面有没有独特的（合胞体来源）群体？
同时检查"融合/黏附"基因（转谷氨酰胺酶、hemocytin、filamin、parvin）在哪个群体表达。
"""
import io, sys, os
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
import scipy.sparse as sp
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"
ORD9 = ["PH", "GH", "MH1", "MH2", "Epi", "F", "R", "B", "NEU"]
BASE = {"2": "PH", "1": "GH", "0": "MH1", "3": "MH2", "4": "Epi", "5": "F", "6": "R", "7": "B", "8": "NE", "9": "Neu"}
# 各组 signature（基因 ID）
four = pd.read_csv(os.path.join(OUT, "blood4_dotplot_genes.csv")); four["id"] = four["gid"].astype(str)
nine = pd.read_csv(os.path.join(OUT, "marker_heatmap_selected_fixed.csv"))
SIG = {g: four[four.grp == g].id.tolist() for g in ["PH", "GH", "MH1", "MH2"]}
for g in ["Epi", "F", "R", "B", "NEU"]:
    SIG[g] = nine[nine.grp == g]["names"].astype(str).tolist()
PANEL = {"TGM(血细胞转谷氨酰胺酶)": ["ncbi_113823945", "ncbi_113801263"],
         "TGM1": ["ncbi_113802929"], "hemocytin/SSPO": ["ncbi_113824171"],
         "ccnd2": ["ncbi_113814652"], "Grik2": ["ncbi_113823802"], "PPAF2": ["ncbi_113812814"]}

raw = ad.read_h5ad(os.path.join(OUT, "adata_raw.h5ad"))
sx = ad.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"), backed="r")
kidx = set(map(str, sx.obs_names))
kept = np.array([str(x) in kidx for x in raw.obs_names])
D = raw[~kept].copy()                      # 10,137 被删细胞
print("被删细胞:", D.shape)
D.X = D.layers["counts"].copy()
sc.pp.filter_genes(D, min_cells=10)
sc.pp.normalize_total(D, target_sum=1e4); sc.pp.log1p(D)
sc.pp.highly_variable_genes(D, n_top_genes=2000); D = D[:, D.var.highly_variable].copy()
sc.pp.scale(D, max_value=10); sc.tl.pca(D, n_comps=30, svd_solver="arpack")
sc.pp.neighbors(D, n_neighbors=15); sc.tl.leiden(D, resolution=0.6, flavor="igraph", n_iterations=2,
                                                 directed=False, key_added="dl")
print("Leiden 群数:", D.obs["dl"].nunique())
# 用保留细胞的 9 群 signature 给每个 Leiden 群打分
K = raw[kept].copy(); K.X = K.layers["counts"].copy()
sc.pp.normalize_total(K, target_sum=1e4); sc.pp.log1p(K)
gpos = {g: i for i, g in enumerate(K.var_names.astype(str))}
sigK = {}
for g, gl in SIG.items():
    gi = [gpos[x] for x in gl if x in gpos]
    sigK[g] = np.asarray(K[:, gi].X.mean(1)).ravel() if gi else np.zeros(K.n_obs)
ct9 = sx.obs["vendor_cluster"].astype(str).map(BASE).replace({"NE": "NEU", "Neu": "NEU"}).reindex(list(map(str, K.obs_names))).values
refmean = pd.DataFrame(sigK).groupby(ct9).mean().reindex(ORD9)
print("\n保留细胞的 9 群 signature 均值（参照）："); print(refmean.round(2).to_string())

D2 = raw[~kept].copy(); D2.X = D2.layers["counts"].copy()
sc.pp.normalize_total(D2, target_sum=1e4); sc.pp.log1p(D2)
gp2 = {g: i for i, g in enumerate(D2.var_names.astype(str))}
S = {}
for g, gl in SIG.items():
    gi = [gp2[x] for x in gl if x in gp2]
    S[g] = np.asarray(D2[:, gi].X.mean(1)).ravel() if gi else np.zeros(D2.n_obs)
for name, gl in PANEL.items():
    gi = [gp2[x] for x in gl if x in gp2]
    S[name] = np.asarray(D2[:, gi].X.mean(1)).ravel() if gi else np.zeros(D2.n_obs)
S["pct_mt"] = D2.obs.pct_counts_mt.values
S["n_genes"] = D2.obs.n_genes_by_counts.values
S["umi"] = D2.obs.total_counts.values
tab = pd.DataFrame(S, index=D2.obs_names); tab["dl"] = D.obs["dl"].values
prof = tab.groupby("dl").mean(numeric_only=True)
prof["n"] = tab.groupby("dl").size()
print("\n=== 被删细胞的 Leiden 群画像（与保留细胞 9 群 signature 对照）===")
print(prof[["n", "umi", "n_genes", "pct_mt"] + ORD9 + list(PANEL)].round(2).to_string())
prof.to_csv(os.path.join(OUT, "dropped_cells_cluster_profile.csv"))
