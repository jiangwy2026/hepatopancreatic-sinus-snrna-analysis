# -*- coding: utf-8 -*-
"""
125_giant_cells_whoarethey.py
(a) 反推本项目 QC 的实际阈值 (b) 巨型条码是环境 RNA 团块还是真细胞 (c) 与 soup 谱的相关性
"""
import io, sys, os
import numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
OUT = r"D:\桌面\single cell\analysis_output"
sym = pd.read_csv(os.path.join(r"D:\桌面\single cell", "ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result",
      "ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result", "1.Expression", "expressions", "blank", "features.tsv.gz"),
      sep="\t", header=None, names=["gid", "sym", "type"])["sym"].astype(str).values

raw = ad.read_h5ad(os.path.join(OUT, "adata_raw.h5ad"))
sx = ad.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"), backed="r")
kept = np.array([str(x) in set(map(str, sx.obs_names)) for x in raw.obs_names])
o = raw.obs
print("=== (a) 本项目保留细胞的实际 QC 区间 ===")
k = o[kept]
print("保留: n_genes %d–%d (中位 %.0f) | pct_mt %.2f–%.2f (中位 %.2f) | total %.0f–%.0f" % (
    k.n_genes_by_counts.min(), k.n_genes_by_counts.max(), k.n_genes_by_counts.median(),
    k.pct_counts_mt.min(), k.pct_counts_mt.max(), k.pct_counts_mt.median(),
    k.total_counts.min(), k.total_counts.max()))
d = o[~kept]
print("丢弃: n_genes %d–%d (中位 %.0f) | pct_mt %.2f–%.2f (中位 %.2f) | total %.0f–%.0f" % (
    d.n_genes_by_counts.min(), d.n_genes_by_counts.max(), d.n_genes_by_counts.median(),
    d.pct_counts_mt.min(), d.pct_counts_mt.max(), d.pct_counts_mt.median(),
    d.total_counts.min(), d.total_counts.max()))
for thr in (8, 10, 12, 15):
    m = (d.pct_counts_mt > thr).sum(); print("   丢弃细胞中 pct_mt > %2d%% : %5d (%.0f%%)" % (thr, m, 100*m/len(d)))

X = raw.layers["counts"].tocsr()
tot = np.asarray(X.sum(1)).ravel()
print("\n=== (b)(c) 巨型条码的成分分析 ===")
# soup 谱（SoupX 反推的环境 RNA 比例），仅 treated
soup_t = np.load(os.path.join(OUT, "_tmp_soup_treated.npy"))
soup_t = soup_t / soup_t.sum()
# PH 群的表达谱（用保留细胞里的 PH 平均）
lab = sx.obs["vendor_cluster"].astype(str).map({"2":"PH","1":"GH","0":"MH1","3":"MH2","4":"Epi","5":"F","6":"R","7":"B","8":"NE","9":"Neu"})
ph_cells = set(lab[lab == "PH"].index)
ph_idx = np.array([str(x) in ph_cells for x in raw.obs_names])
Xp = X[ph_idx]; php = np.asarray(Xp.sum(0)).ravel().astype(float); php /= php.sum()

for name, i in [("treated_115k", int(np.where(np.array(raw.obs_names) == "treated_" + str(o.index[(raw.obs_names == "treated")]))[0][0]) if False else None)]:
    pass
# 直接取 UMI 最大的两个 + 若干 PH 高 UMI
order = np.argsort(-tot)
top = order[:6]
rows = []
for i in list(top) + [int(np.where(raw.obs_names == "treated_CAACATTAGGTAATCT")[0][0]) if "treated_CAACATTAGGTAATCT" in list(raw.obs_names) else -1]:
    if i < 0: continue
    v = np.asarray(X[i].todense()).ravel().astype(float)
    if v.sum() == 0: continue
    p = v / v.sum()
    r_soup = float(np.corrcoef(p, soup_t)[0, 1])
    r_ph = float(np.corrcoef(p, php)[0, 1])
    topg = np.argsort(-p)[:8]
    rows.append({"cell": str(raw.obs_names[i]), "umi": int(tot[i]), "n_genes": int((v > 0).sum()),
                 "pct_mt": round(float(o.pct_counts_mt.values[i]), 2),
                 "r_soup": round(r_soup, 3), "r_PH": round(r_ph, 3),
                 "top_genes": ", ".join("%s(%.1f%%)" % (sym[g], 100 * p[g]) for g in topg[:5])})
print(pd.DataFrame(rows).to_string(index=False))
# 参照：保留细胞中普通 PH 细胞的这两个相关性
sub = np.where(ph_idx & kept)[0][:200]
rs, rp = [], []
for i in sub:
    v = np.asarray(X[i].todense()).ravel().astype(float); p = v / v.sum()
    rs.append(np.corrcoef(p, soup_t)[0, 1]); rp.append(np.corrcoef(p, php)[0, 1])
print("\n参照（200 个 PH 细胞）: r_soup 中位 %.3f | r_PH 中位 %.3f" % (np.median(rs), np.median(rp)))
