# -*- coding: utf-8 -*-
"""
132_pertype_qc.py (step 1)
把 19,795 个 CellRanger 细胞全部标上类型（kNN 标签转移），再按类型定 QC 阈值。
"""
import io, sys, os
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
import scipy.sparse as sp
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"
root = r"D:\桌面\single cell"
ORD = ["PH", "GH", "MH1", "MH2", "Epi", "F", "R", "B", "NEU"]
BASE = {"2": "PH", "1": "GH", "0": "MH1", "3": "MH2", "4": "Epi", "5": "F", "6": "R", "7": "B", "8": "NEU", "9": "NEU"}

cl = pd.read_csv(os.path.join(root, "ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result",
    "ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result", "2.Cluster", "2.cluster", "Cells.cluster.list.xls"), sep="\t")
cl["ct"] = cl["Cluster"].astype(str).map(BASE)
lab = dict(zip(cl["Cells"].astype(str), cl["ct"]))

A = ad.read_h5ad(os.path.join(OUT, "adata_raw.h5ad"))      # 19,795
A.X = A.layers["counts"].copy()
lab0 = np.array([lab.get(str(x), None) for x in A.obs_names], dtype=object)
print("已有标签: %d | 待转移: %d" % ((lab0 != None).sum(), (lab0 == None).sum()))

# kNN 标签转移：先用全部细胞做 PCA（用 HVG）
sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
sc.pp.highly_variable_genes(A, n_top_genes=2000)
Ah = A[:, A.var.highly_variable].copy()
sc.pp.scale(Ah, max_value=10); sc.tl.pca(Ah, n_comps=30, svd_solver="arpack")
sc.pp.neighbors(Ah, n_neighbors=15, use_rep="X_pca")
conn = Ah.obsp["connectivities"]
known = np.where(lab0 != None)[0]; unknown = np.where(lab0 == None)[0]
K = conn[unknown][:, known]                      # 未知 x 已知 的连接
W = np.asarray(K.todense()) if sp.issparse(K) else np.asarray(K)
print("标签转移: 每个未知细胞的邻居数中位 %.0f" % np.median(W.sum(1)))
kidx = {g: i for i, g in enumerate(ORD)}
lidx = np.array([kidx[str(lab0[i])] for i in known])
votes = np.zeros((len(unknown), len(ORD)))
for j in range(len(ORD)):
    votes[:, j] = W[:, lidx == j].sum(1)
pred = np.array(ORD)[votes.argmax(1)]
conf = votes.max(1) / np.maximum(votes.sum(1), 1e-9)
lab_all = lab0.copy(); lab_all[unknown] = pred
print("转移后各组细胞数:")
print(pd.Series(lab_all).value_counts().reindex(ORD).to_string())
print("转移置信度(最高票占比): 中位 %.2f, <25%% 的细胞 %d 个" % (np.median(conf), (conf < 0.25).sum()))

obs = A.obs.copy(); obs["ct"] = lab_all
obs.to_csv(os.path.join(OUT, "all_19795_celltype.csv"))

print("\n=== 全 19,795 细胞按类型分组的 QC 分布 ===")
rows = []
for g in ORD:
    d = obs[obs["ct"] == g]
    mt = d.pct_counts_mt.values; ng = d.n_genes_by_counts.values
    med, mad = np.median(mt), np.median(np.abs(mt - np.median(mt))) * 1.4826
    rows.append({"ct": g, "n": len(d), "mt_med": round(med, 2), "mt_mad": round(mad, 2),
                 "mt_p90": round(np.percentile(mt, 90), 1), "mt_p95": round(np.percentile(mt, 95), 1),
                 "mt_cut3": round(med + 3 * mad, 1), "mt_cut5": round(med + 5 * mad, 1),
                 "ng_med": int(np.median(ng)), "keep_3mad": int(((mt < med + 3 * mad) & (ng >= 200)).sum()),
                 "keep_5mad": int(((mt < med + 5 * mad) & (ng >= 200)).sum()),
                 "keep_p99": int(((mt < np.percentile(mt, 99)) & (ng >= 200)).sum())})
T = pd.DataFrame(rows)
print(T.to_string(index=False))
T.to_csv(os.path.join(OUT, "pertype_qc_all19795.csv"), index=False)
print("\n合计保留: 3MAD规则 %d | 5MAD规则 %d | p99规则 %d （厂商为 9,658）" % (
    T.keep_3mad.sum(), T.keep_5mad.sum(), T.keep_p99.sum()))
print("\n按样本×类型统计（5MAD 规则）:")
keep5 = np.zeros(len(obs), bool)
for g in ORD:
    m = (obs["ct"] == g).values
    mt = obs.pct_counts_mt.values; ng = obs.n_genes_by_counts.values
    med, mad = np.median(mt[m]), np.median(np.abs(mt[m] - np.median(mt[m]))) * 1.4826
    keep5 |= (m & (mt < med + 5 * mad) & (ng >= 200))
obs2 = obs.assign(keep5=keep5)
print(pd.crosstab(obs2["ct"], obs2["sample"], values=obs2["keep5"], aggfunc="sum").reindex(ORD).to_string())
print("\n各样本总保留:", obs2.groupby("sample")["keep5"].sum().to_dict())
