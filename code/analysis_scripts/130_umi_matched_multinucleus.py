# -*- coding: utf-8 -*-
"""
130_umi_matched_multinucleus.py
S/G2M 双阳性比例在 0h vs 6h 的差异，是不是因为 0h 测到更多 RNA（深度效应）？
做 UMI 分层匹配后重新比较；并给出"多核样"的严格定义。
"""
import io, sys, os
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
OUT = r"D:\桌面\single cell\analysis_output"
four = pd.read_csv(os.path.join(OUT, "blood4_dotplot_genes.csv")); four["id"] = four["gid"].astype(str)
SIGPH = four[four.grp == "PH"].id.tolist()
S_GENES = ["PCNA","MCM6","CDC6","ORC1","ORC4","ORC5","POLA2","PRIM1","FEN1","TYMS","DHFR","HIRA","CHAF1A"]
G2M_GENES = ["CCNB3","AURKA","TOP2A","TPX2","ASPM","Pka-C1","KIF20A","KIF23","RACGAP1","NDC80","INCENP","PRC1","Patronin"]

raw = ad.read_h5ad(os.path.join(OUT, "adata_raw.h5ad"))
X = raw.layers["counts"].tocsr()
tot = np.asarray(X.sum(1)).ravel()
norm = sp.diags(1e4 / np.maximum(tot, 1)) @ X
ln = norm.copy(); ln.data = np.log1p(ln.data)
ids = raw.var_names.astype(str).values
gp = {g: i for i, g in enumerate(ids)}
sym = raw.var["symbol"].astype(str).values
tg = np.asarray(X.sum(0)).ravel()
best = {}
for i, s in enumerate(sym):
    if s not in best or tg[i] > tg[best[s]]: best[s] = i
def mid(gl): return [gp[x] for x in gl if x in gp]
def msym(gl): return [gp[str(ids[best[x]])] for x in gl if x in best]
def mexp(gi): return np.asarray(ln[:, gi].mean(1)).ravel()
T = pd.DataFrame({
    "sample": raw.obs["sample"].astype(str).values,
    "umi": tot, "n_genes": np.asarray((X > 0).sum(1)).ravel(),
    "pct_mt": raw.obs.pct_counts_mt.values,
    "PH_sig": mexp(mid(SIGPH)), "S": mexp(msym(S_GENES)), "G2M": mexp(msym(G2M_GENES)),
})
T["double_pos"] = (T.S > 0) & (T.G2M > 0)
q = T.PH_sig.quantile(0.90)
P = T[T.PH_sig >= q].copy()
print("PH 样细胞（PH_sig ≥ %.2f）: blank %d | treated %d" % (q, (P["sample"]=="blank").sum(), (P["sample"]=="treated").sum()))
print("\n=== 未匹配：双阳性比例 ===")
print(P.groupby("sample")["double_pos"].agg(["mean", "size"]).assign(pct=lambda d: (100*d["mean"]).round(1)).to_string())

print("\n=== UMI 分层匹配后 ===")
P["umbin"] = pd.cut(P.umi, [0, 600, 900, 1300, 2000, 1e9], labels=["<600", "600-900", "900-1300", "1300-2000", ">2000"])
tab = P.groupby(["umbin", "sample"], observed=True)["double_pos"].agg(["mean", "size"]).reset_index()
tab["pct"] = (100 * tab["mean"]).round(1)
print(tab.pivot(index="umbin", columns="sample", values=["pct", "size"]).to_string())
print("\n带 S/G2M 分数值的比较（仅 UMI 匹配的分层）:")
print(P.groupby(["umbin", "sample"], observed=True)[["S", "G2M"]].mean().round(3).to_string())

print("\n=== 严格多核样定义 ===")
print("S>0.05 且 G2M>0.15（两者都实质表达）:")
strict = (P.S > 0.05) & (P.G2M > 0.15)
for s in ("blank", "treated"):
    d = P[P["sample"] == s]
    print("   %-8s %5d/%5d = %.2f%%" % (s, strict[d.index].sum(), len(d), 100 * strict[d.index].mean()))
sub = P[P.umbin.isin(["900-1300", "1300-2000"])]
st = (sub.S > 0.05) & (sub.G2M > 0.15)
for s in ("blank", "treated"):
    d = sub[sub["sample"] == s]
    print("   [UMI 900–2000 匹配] %-8s %4d/%4d = %.2f%%" % (s, st[d.index].sum(), len(d), 100 * st[d.index].mean()))
