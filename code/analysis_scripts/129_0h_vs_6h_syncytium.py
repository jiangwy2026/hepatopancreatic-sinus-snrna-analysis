# -*- coding: utf-8 -*-
"""
129_0h_vs_6h_syncytium.py
blank = 0 h（合胞体完整） ; treated = 6 h（已 cellularization）
如果合胞体被捕获：应在 0 h 富集"大型/多核样"细胞。逐项检验。
"""
import io, sys, os
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
import scipy.sparse as sp
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"
four = pd.read_csv(os.path.join(OUT, "blood4_dotplot_genes.csv")); four["id"] = four["gid"].astype(str)
SIGPH = four[four.grp == "PH"].id.tolist()
PANEL = {"TGM": ["ncbi_113823945", "ncbi_113801263"], "hemocytin": ["ncbi_113824171"],
         "ccnd2": ["ncbi_113814652"], "Grik2": ["ncbi_113823802"], "PPAF2": ["ncbi_113812814"],
         "Igfbp7": ["ncbi_113811702"], "SSPO": ["ncbi_113824171"]}
S_GENES = ["PCNA","MCM6","CDC6","ORC1","ORC4","ORC5","POLA2","PRIM1","FEN1","TYMS","DHFR","HIRA","CHAF1A"]
G2M_GENES = ["CCNB3","AURKA","TOP2A","TPX2","ASPM","Pka-C1","KIF20A","KIF23","RACGAP1","NDC80","INCENP","PRC1","Patronin"]

raw = ad.read_h5ad(os.path.join(OUT, "adata_raw.h5ad"))
X = raw.layers["counts"].tocsr()
tot = np.asarray(X.sum(1)).ravel(); ng = np.asarray((X > 0).sum(1)).ravel()
sample = raw.obs["sample"].astype(str).values
norm = sp.diags(1e4 / np.maximum(tot, 1)) @ X
ln = norm.copy(); ln.data = np.log1p(ln.data)
var_ids = raw.var_names.astype(str).values
gp = {g: i for i, g in enumerate(var_ids)}
sym = (raw.var["symbol"] if "symbol" in raw.var.columns else raw.var["sym"]).astype(str).values
tot_gene = np.asarray(X.sum(0)).ravel()
best = {}
for i, s in enumerate(sym):
    if s not in best or tot_gene[i] > tot_gene[best[s]]: best[s] = i
def ids(gl): return [gp[str(var_ids[best[x]])] for x in gl if x in best]      # 输入基因名
def ids_by_id(gl): return [gp[x] for x in gl if x in gp]                       # 输入基因 ID
def mean_expr(gl):
    gi = ids(gl)
    return np.asarray(ln[:, gi].mean(1)).ravel() if gi else np.zeros(len(sample))
def mean_expr_id(gl):
    gi = ids_by_id(gl)
    return np.asarray(ln[:, gi].mean(1)).ravel() if gi else np.zeros(len(sample))

S = {}
S["PH_sig"] = mean_expr_id(SIGPH)
for k, v in PANEL.items(): S[k] = mean_expr_id(v)
S["S_score"] = mean_expr(S_GENES); S["G2M_score"] = mean_expr(G2M_GENES)
T = pd.DataFrame(S); T["sample"] = sample; T["umi"] = tot; T["n_genes"] = ng
T["pct_mt"] = raw.obs.pct_counts_mt.values

print("=== (a) 全 19,795 细胞：0h vs 6h 的 RNA 含量分布 ===")
for s in ("blank", "treated"):
    v = T.umi[T["sample"] == s]
    print("  %-8s n=%5d | UMI 中位 %.0f p95 %.0f p99 %.0f 最大 %d | >3000: %d (%.2f%%) | >6000: %d (%.2f%%)" % (
        s, len(v), np.median(v), np.percentile(v, 95), np.percentile(v, 99), v.max(),
        (v > 3000).sum(), 100 * (v > 3000).mean(), (v > 6000).sum(), 100 * (v > 6000).mean()))
print("  %-8s n_genes 中位 %.0f | >1500: %d (%.2f%%)" % ("blank", T[T["sample"]=="blank"].n_genes.median(),
      (T[T["sample"]=="blank"].n_genes > 1500).sum(), 100*(T[T["sample"]=="blank"].n_genes > 1500).mean()))
print("  %-8s n_genes 中位 %.0f | >1500: %d (%.2f%%)" % ("treated", T[T["sample"]=="treated"].n_genes.median(),
      (T[T["sample"]=="treated"].n_genes > 1500).sum(), 100*(T[T["sample"]=="treated"].n_genes > 1500).mean()))

print("\n=== (b) PH 样细胞（PH signature 前 10%）在 0h vs 6h 的对比 ===")
q = T.PH_sig.quantile(0.90)
ph = T[T.PH_sig >= q].copy()
print("  阈值 PH_sig >= %.2f，共 %d 个细胞" % (q, len(ph)))
for s in ("blank", "treated"):
    d = ph[ph["sample"] == s]
    print("  %-8s n=%4d | UMI 中位 %.0f p95 %.0f | n_genes 中位 %.0f | pct_mt 中位 %.1f" % (
        s, len(d), d.umi.median(), np.percentile(d.umi, 95), d.n_genes.median(), d.pct_mt.median()))
    print("            TGM %.2f | hemocytin %.2f | ccnd2 %.2f | Grik2 %.2f | PPAF2 %.2f | S %.3f | G2M %.3f" % (
        d.TGM.mean(), d.hemocytin.mean(), d.ccnd2.mean(), d.Grik2.mean(), d.PPAF2.mean(), d.S_score.mean(), d.G2M_score.mean()))
    print("            S>0 且 G2M>0 双阳性比例 %.1f%% | S>0 %.1f%% | G2M>0 %.1f%%" % (
        100 * ((d.S_score > 0) & (d.G2M_score > 0)).mean(), 100 * (d.S_score > 0).mean(), 100 * (d.G2M_score > 0).mean()))

print("\n=== (c) '大型/多核样'候选（UMI > 3000 且 PH_sig 高）在两组中的分布 ===")
cand = T[(T.umi > 3000) & (T.PH_sig >= q)]
print("  合计 %d 个 | blank %d (%.2f%% of blank) | treated %d (%.2f%% of treated)" % (
    len(cand), (cand["sample"] == "blank").sum(), 100 * (cand["sample"] == "blank").mean() * len(cand) / (sample == "blank").sum(),
    (cand["sample"] == "treated").sum(), 100 * (cand["sample"] == "treated").mean() * len(cand) / (sample == "treated").sum()))
print(cand.sort_values("umi", ascending=False).head(12)[["sample","umi","n_genes","pct_mt","PH_sig","TGM","G2M_score"]].round(2).to_string())
