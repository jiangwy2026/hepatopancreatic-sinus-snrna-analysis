# -*- coding: utf-8 -*-
"""131_pertype_qc_explore.py - 各细胞类型的 QC 分布（用于按类型定阈值）"""
import io, sys, os
import numpy as np, pandas as pd, anndata as ad
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
OUT = r"D:\桌面\single cell\analysis_output"
root = r"D:\桌面\single cell"
cl = pd.read_csv(os.path.join(root, "ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result",
    "ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result", "2.Cluster", "2.cluster", "Cells.cluster.list.xls"),
    sep="\t")
cl["Cells"] = cl["Cells"].astype(str); cl["Cluster"] = cl["Cluster"].astype(str)
print("厂商 cell-cluster 表:", cl.shape, "| 列:", list(cl.columns), "| 样本:", cl["Cells"].str.split("_").str[0].value_counts().to_dict())
BASE = {"2": "PH", "1": "GH", "0": "MH1", "3": "MH2", "4": "Epi", "5": "F", "6": "R", "7": "B", "8": "NEU", "9": "NEU"}
cl["ct"] = cl["Cluster"].map(BASE)
lab = cl.set_index("Cells")["ct"]
A = ad.read_h5ad(os.path.join(OUT, "adata_raw.h5ad"))
A.obs["ct"] = lab.reindex(list(map(str, A.obs_names))).values
print("未能标到类型的细胞:", int(A.obs["ct"].isna().sum()))
print("\n各组细胞数:\n", A.obs["ct"].value_counts(dropna=False).to_string())
ORD = ["PH", "GH", "MH1", "MH2", "Epi", "F", "R", "B", "NEU"]
rows = []
for g in ORD:
    d = A.obs[A.obs["ct"] == g]
    row = {"ct": g, "n": len(d)}
    for col, tag in (("pct_counts_mt", "mt"), ("n_genes_by_counts", "ng"), ("total_counts", "umi")):
        v = d[col].values
        med = np.median(v); mad = np.median(np.abs(v - med)) * 1.4826
        row.update({tag + "_med": round(med, 2), tag + "_mad": round(mad, 2),
                    tag + "_p90": round(np.percentile(v, 90), 2), tag + "_p95": round(np.percentile(v, 95), 2),
                    tag + "_p99": round(np.percentile(v, 99), 2), tag + "_max": round(v.max(), 1),
                    tag + "_cut3mad": round(med + 3 * mad, 2), tag + "_cut5mad": round(med + 5 * mad, 2)})
    rows.append(row)
t = pd.DataFrame(rows)
print("\n=== 各类型 QC 分布 ===")
print(t[["ct", "n", "mt_med", "mt_mad", "mt_p90", "mt_p95", "mt_p99", "mt_max", "mt_cut3mad", "mt_cut5mad"]].to_string(index=False))
print()
print(t[["ct", "ng_med", "ng_p10", "um = ng_p90".split(" = ")[0] if False else "ng_p90", "ng_p99"]].to_string(index=False) if False else t[["ct", "ng_med", "ng_p90", "ng_p99", "umi_med", "umi_p99"]].to_string(index=False))
# 每类型：用不同规则保留多少
print("\n=== 不同阈值规则下各类型保留比例（% 全部细胞）===")
out = []
for g in ORD:
    d = A.obs[A.obs["ct"] == g]
    n = len(d)
    mt = d.pct_counts_mt.values; ng = d.n_genes_by_counts.values
    med = np.median(mt); mad = np.median(np.abs(mt - med)) * 1.4826
    r = {"ct": g, "n": n,
         "全局 mt<10 & ng>=200": 100 * ((mt < 10) & (ng >= 200)).mean(),
         "mt<med+3MAD & ng>=200": 100 * ((mt < med + 3 * mad) & (ng >= 200)).mean(),
         "mt<med+5MAD & ng>=200": 100 * ((mt < med + 5 * mad) & (ng >= 200)).mean(),
         "mt<p95 & ng>=200": 100 * ((mt < np.percentile(mt, 95)) & (ng >= 200)).mean(),
         "仅 ng>=200": 100 * (ng >= 200).mean()}
    out.append(r)
print(pd.DataFrame(out).round(1).to_string(index=False))
t.to_csv(os.path.join(OUT, "pertype_qc_stats.csv"), index=False)
