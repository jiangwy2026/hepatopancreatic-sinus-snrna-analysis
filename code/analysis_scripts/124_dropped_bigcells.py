# -*- coding: utf-8 -*-
"""
124_dropped_bigcells.py
被丢弃的高 UMI 细胞到底是什么？是死细胞、环境 RNA 团块、双细胞，还是"合胞体/合胞体解离产物"？
"""
import io, sys, os
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
import scipy.sparse as sp
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"
ORD = ["PH", "GH", "MH1", "MH2", "Epi", "F", "R", "B", "NEU"]
BASE = {"2": "PH", "1": "GH", "0": "MH1", "3": "MH2", "4": "Epi", "5": "F", "6": "R", "7": "B", "8": "NE", "9": "Neu"}
# 各组 marker（用基因 ID，来自已修正的 marker 表）
four = pd.read_csv(os.path.join(OUT, "blood4_dotplot_genes.csv")); four["id"] = four["gid"].astype(str)
nine = pd.read_csv(os.path.join(OUT, "marker_heatmap_selected_fixed.csv"))
SIG = {g: four[four.grp == g].id.tolist() for g in ["PH", "GH", "MH1", "MH2"]}
for g in ["Epi", "F", "R", "B", "NEU"]:
    d = nine[nine.grp == g]
    SIG[g] = d["names"].astype(str).tolist() if "names" in d.columns else []
print("各组 signature 基因数:", {k: len(v) for k, v in SIG.items()})

raw = ad.read_h5ad(os.path.join(OUT, "adata_raw.h5ad"))
sx = ad.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"), backed="r")
kept = set(map(str, sx.obs_names))
raw.obs["kept"] = [str(x) in kept for x in raw.obs_names]
X = raw.layers["counts"].copy()
X = X / np.maximum(np.asarray(X.sum(1)), 1) * 1e4
X = sp.csr_matrix(np.log1p(X))
gpos = {g: i for i, g in enumerate(raw.var_names.astype(str))}

def score(idx_cells):
    s = {}
    for g, gl in SIG.items():
        gi = [gpos[x] for x in gl if x in gpos]
        if not gi: continue
        s[g] = np.asarray(X[idx_cells][:, gi].mean(1)).ravel()
    return pd.DataFrame(s, index=np.array(raw.obs_names)[idx_cells])

# 参照：本项目保留细胞里各群的 signature 均值（作为"正常"基线）
kept_idx = np.where(raw.obs["kept"].values)[0]
lab = sx.obs["vendor_cluster"].astype(str).map(BASE).replace({"NE": "NEU", "Neu": "NEU"})
lab = lab.reindex([str(x) for x in raw.obs_names[kept_idx]]).values
ref = score(kept_idx); ref["ct"] = lab
base = ref.groupby("ct")[ORD].mean()
print("\n=== 本项目保留细胞中各群的 signature 平均得分（基线）===")
print(base.round(3).to_string())

big = np.where(raw.obs.total_counts.values > 10000)[0]
print("\n=== UMI > 10,000 的被丢弃细胞（%d 个）===" % len(big))
sb = score(big)
o = raw.obs.iloc[big]
out = pd.DataFrame({
    "sample": o["sample"].values, "umi": o["total_counts"].values, "n_genes": o["n_genes_by_counts"].values,
    "pct_mt": o["pct_counts_mt"].round(1).values,
    "best": sb[ORD].idxmax(1).values, "best_score": sb[ORD].max(1).round(2).values,
})
print(out.to_string(index=False))
print("\n这些细胞的 signature 得分矩阵：")
print(sb.round(2).to_string())

mid = np.where((raw.obs.total_counts.values > 5000) & (raw.obs.total_counts.values <= 10000) & (~raw.obs["kept"].values))[0]
sm = score(mid)
print("\n=== UMI 5,000–10,000 被丢弃细胞（%d 个）汇总 ===" % len(mid))
print(pd.DataFrame({"best": sm[ORD].idxmax(1).values}).value_counts().to_string())
print("平均 pct_mt: %.1f%% | 平均 n_genes: %.0f" % (
    raw.obs.iloc[mid].pct_counts_mt.mean(), raw.obs.iloc[mid].n_genes_by_counts.mean()))
json_dump = {"n_big": int(len(big)), "n_mid": int(len(mid)),
             "big_best_counts": pd.Series(sb[ORD].idxmax(1).values).value_counts().to_dict()}
import json; json.dump(json_dump, open(os.path.join(OUT, "DROPPED_BIGCELLS.json"), "w"), indent=1)
