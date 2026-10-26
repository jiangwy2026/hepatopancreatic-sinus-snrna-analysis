# -*- coding: utf-8 -*-
"""
127_syncytium_tests.py
(A) 被删的 PH 样群(dl 11) vs 保留 PH：祖细胞程序、细胞周期、样本构成
(B) 环境 RNA（soup）中是否富集 PH/祖细胞程序基因 -> 检验"合胞体裂解成环境 RNA"假说
"""
import io, sys, os
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
import scipy.sparse as sp
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"
ORD9 = ["PH", "GH", "MH1", "MH2", "Epi", "F", "R", "B", "NEU"]
BASE = {"2": "PH", "1": "GH", "0": "MH1", "3": "MH2", "4": "Epi", "5": "F", "6": "R", "7": "B", "8": "NE", "9": "Neu"}
four = pd.read_csv(os.path.join(OUT, "blood4_dotplot_genes.csv")); four["id"] = four["gid"].astype(str)
nine = pd.read_csv(os.path.join(OUT, "marker_heatmap_selected_fixed.csv"))
SIG = {g: four[four.grp == g].id.tolist() for g in ["PH", "GH", "MH1", "MH2"]}
for g in ["Epi", "F", "R", "B", "NEU"]:
    SIG[g] = nine[nine.grp == g]["names"].astype(str).tolist()
PANEL = {"TGM": ["ncbi_113823945", "ncbi_113801263"], "hemocytin": ["ncbi_113824171"],
         "ccnd2": ["ncbi_113814652"], "Grik2": ["ncbi_113823802"], "PPAF2": ["ncbi_113812814"]}
S_GENES = ["PCNA","MCM6","CDC6","ORC1","ORC4","ORC5","POLA2","PRIM1","FEN1","TYMS","DHFR","HIRA","CHAF1A"]
G2M_GENES = ["CCNB3","AURKA","TOP2A","TPX2","ASPM","Pka-C1","KIF20A","KIF23","RACGAP1","NDC80","INCENP","PRC1","Patronin"]

raw = ad.read_h5ad(os.path.join(OUT, "adata_raw.h5ad"))
sx = ad.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"), backed="r")
kidx = set(map(str, sx.obs_names))
kept = np.array([str(x) in kidx for x in raw.obs_names])
ct9 = sx.obs["vendor_cluster"].astype(str).map(BASE).replace({"NE": "NEU", "Neu": "NEU"}).reindex(list(map(str, raw.obs_names[kept]))).values

# ---------- (A) 保留细胞 + 被删细胞的同一套指标
drops = pd.read_csv(os.path.join(OUT, "dropped_cells_cluster_profile.csv"), index_col=0)
dl = pd.read_pickle(os.path.join(OUT, "_tmp_soup_store.pkl")) if False else None
X = raw.layers["counts"].tocsr()
tot = np.asarray(X.sum(1)).ravel()
norm = sp.diags(1e4 / np.maximum(tot, 1)) @ X
ln = norm.copy(); ln.data = np.log1p(ln.data)
gp = {g: i for i, g in enumerate(raw.var_names.astype(str))}
# 符号 -> 取总 counts 最高的那个基因 ID（避免同名基因取错）
_sym = (raw.var["symbol"] if "symbol" in raw.var.columns else raw.var["sym"]).astype(str).values
_tot_gene = np.asarray(X.sum(0)).ravel()
_best = {}
for i, s in enumerate(_sym):
    if s not in _best or _tot_gene[i] > _tot_gene[_best[s]]:
        _best[s] = i
def _ids(gl):
    return [gp[str(raw.var_names[_best[x]])] for x in gl if x in _best]
def panel_scores(idx):
    out = {}
    for g, gl in SIG.items():
        gi = [gp[x] for x in gl if x in gp]
        out[g] = np.asarray(ln[idx][:, gi].mean(1)).ravel()
    for name, gl in PANEL.items():
        gi = [gp[x] for x in gl if x in gp]
        out[name] = np.asarray(ln[idx][:, gi].mean(1)).ravel()
    for nm, gl in (("S", S_GENES), ("G2M", G2M_GENES)):
        gi = _ids(gl)
        out[nm] = np.asarray(ln[idx][:, gi].mean(1)).ravel()
    return pd.DataFrame(out, index=np.array(raw.obs_names)[idx])
ki = np.where(kept)[0]
K = panel_scores(ki); K["ct"] = ct9; K["sample"] = raw.obs["sample"].values[ki]
print("=== 保留细胞各群指标（含祖细胞程序与细胞周期）===")
print(K.groupby("ct")[["PH","TGM","hemocytin","ccnd2","Grik2","PPAF2","S","G2M"]].mean().reindex(ORD9).round(2).to_string())
print("\n保留 PH：细胞周期双阳性（S>0 且 G2M>0）比例 %.1f%%" % (
    100 * ((K[K.ct=="PH"].S > 0) & (K[K.ct=="PH"].G2M > 0)).mean()))

# 被删细胞中的 PH 样群 = 126 脚本里的 dl 11（n=332）：重建该群索引
D = raw[~kept].copy(); D.X = D.layers["counts"].copy()
sc.pp.filter_genes(D, min_cells=10)
sc.pp.normalize_total(D, target_sum=1e4); sc.pp.log1p(D)
sc.pp.highly_variable_genes(D, n_top_genes=2000); D = D[:, D.var.highly_variable].copy()
sc.pp.scale(D, max_value=10); sc.tl.pca(D, n_comps=30, svd_solver="arpack")
sc.pp.neighbors(D, n_neighbors=15); sc.tl.leiden(D, resolution=0.6, flavor="igraph", n_iterations=2, directed=False, key_added="dl")
di = np.where(~kept)[0]
Dp = panel_scores(di); Dp["dl"] = D.obs["dl"].values; Dp["sample"] = raw.obs["sample"].values[di]
print("\n=== 被删细胞里的 PH 样群（dl 11）与被删其他群对比 ===")
show = Dp.groupby("dl")[["PH","TGM","hemocytin","ccnd2","Grik2","ppaf2" if "ppaf2" in Dp else "PPAF2","S","G2M"]].mean().round(2)
show["n"] = Dp.groupby("dl").size()
print(show.to_string())
dl11 = Dp.index[Dp.dl == "11"]
print("\ndl 11 的样本构成:", pd.Series(raw.obs["sample"].values[di][Dp.dl.values == "11"]).value_counts().to_dict())
print("dl 11 在 blank/treated 中的占该样本细胞数比例: blank %.2f%%  treated %.2f%%" % (
    100 * (Dp["sample"].values[Dp.dl.values == "11"] == "blank").sum() / (raw.obs["sample"].values[di] == "blank").sum(),
    100 * (Dp["sample"].values[Dp.dl.values == "11"] == "treated").sum() / (raw.obs["sample"].values[di] == "treated").sum()))
print("保留 PH 在样本中的比例: blank %.2f%%  treated %.2f%%" % (
    100 * ((K.ct == "PH") & (K["sample"] == "blank")).sum() / (K["sample"] == "blank").sum(),
    100 * ((K.ct == "PH") & (K["sample"] == "treated")).sum() / (K["sample"] == "treated").sum()))

# ---------- (B) 环境 RNA 是否富集 PH/祖细胞程序
print("\n=== (B) 环境 RNA（SoupX soup 谱）里的基因富集 ===")
F = pd.read_csv(os.path.join(r"D:\桌面\single cell", "ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result",
     "ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result", "1.Expression", "expressions", "blank", "features.tsv.gz"),
     sep="\t", header=None, names=["gid","sym","type"])
syms = F["sym"].astype(str).values
for s in ("blank", "treated"):
    soup = np.load(os.path.join(OUT, "_tmp_soup_%s.npy" % s)).astype(float)
    soup = soup / soup.sum()
    idxs = np.where(kept & (raw.obs["sample"].values == s))[0]
    cellbulk = np.asarray(X[idxs].sum(0)).ravel().astype(float); cellbulk = cellbulk / cellbulk.sum()
    ok = cellbulk > 0
    enr = np.full(len(soup), np.nan); enr[ok] = soup[ok] / cellbulk[ok]
    rows = []
    for g, gl in SIG.items():
        gi = [gp[x] for x in gl if x in gp and ok[gp[x]]]
        if gi: rows.append({"signature": g, "median_soup_enrichment": round(float(np.median(enr[gi])), 2)})
    for name, gl in PANEL.items():
        gi = [gp[x] for x in gl if x in gp and ok[gp[x]]]
        if gi: rows.append({"signature": name, "median_soup_enrichment": round(float(np.median(enr[gi])), 2)})
    # 全局参考：所有基因的中位富集
    allmed = float(np.nanmedian(enr[ok]))
    print("\n[%s] 全部基因 soup/cell 富集中位 = %.2f" % (s, allmed))
    t = pd.DataFrame(rows).sort_values("median_soup_enrichment", ascending=False)
    print(t.to_string(index=False))
    # 最富集的前 15 个基因
    top = np.argsort(-np.nan_to_num(enr, nan=0))[:15]
    print("  soup 中最富集的 15 个基因:", ", ".join("%s(%.1fx)" % (syms[i], enr[i]) for i in top if not np.isnan(enr[i])))
