import sys, io, os
import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
from scipy.stats import rankdata
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace"); pd.set_option("display.width", 250)
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"; fig = out + r"\figures_traj"; os.makedirs(fig, exist_ok=True)
A = sc.read_h5ad(out + r"\adata_hemocyte_lineage.h5ad")
sub = A.obs["sub"].astype(str).values
P = A.obsp["connectivities"].tocsr().astype(float)
P = sp.diags(1.0 / np.asarray(P.sum(1)).ravel()) @ P
m2 = (sub == "2"); m0 = (sub == "0"); m1 = (sub == "1")
R, A0, A1 = np.where(m2)[0], np.where(m0)[0], np.where(m1)[0]
print("增殖亚群 R=%d, 终态 A0(亚群0)=%d, A1(亚群1)=%d" % (len(R), len(A0), len(A1)))
PRR = P[R][:, R].toarray(); PR0 = np.asarray(P[R][:, A0].sum(1)).ravel(); PR1 = np.asarray(P[R][:, A1].sum(1)).ravel()
I = np.eye(len(R))
X = np.linalg.solve(I - PRR, np.column_stack([PR0, PR1]))
p0, p1 = X[:, 0], X[:, 1]
print("\n=== 增殖祖细胞的命运概率（均值）===")
print("  → 亚群0（静息/未分化态）: %.3f" % p0.mean())
print("  → 亚群1（成熟效应血细胞）: %.3f" % p1.mean())
A.obs.loc[A.obs_names[R], "fate_to_0"] = p0
A.obs.loc[A.obs_names[R], "fate_to_1"] = p1
A.obs["fate_bias"] = np.nan
A.obs.loc[A.obs_names[R], "fate_bias"] = p1 - p0
print("\n=== 按组的命运概率 ===")
print(A.obs.loc[A.obs_names[R]].groupby("group")[["fate_to_0","fate_to_1"]].mean().round(3).to_string())
print("\n=== 按拟时序分位（祖细胞内部先后）===")
tmp = A.obs.loc[A.obs_names[R], ["dpt_pseudotime","fate_to_0","fate_to_1"]].copy()
tmp["pt_bin"] = pd.qcut(tmp["dpt_pseudotime"], 4, labels=["Q1最早","Q2","Q3","Q4最晚"])
print(tmp.groupby("pt_bin", observed=True)[["fate_to_0","fate_to_1"]].mean().round(3).to_string())
# 命运偏向的驱动基因
Xg = A.X if not hasattr(A.X,"todense") else np.asarray(A.X.todense())
Xr = Xg[R]
bias = p1 - p0; rk = rankdata(bias); rk = (rk - rk.mean())/(rk.std()+1e-12)
var = Xr.var(0); hv = np.argsort(-var)[:4000]
rx = rankdata(Xr[:, hv], axis=0); rx = (rx - rx.mean(0))/(rx.std(0)+1e-12)
rho = (rx * rk[:, None]).mean(0)
res = pd.DataFrame({"gene": A.var["symbol"].values[hv].astype(str), "rho": rho}).sort_values("rho")
res.to_csv(out + r"\fate_bias_genes.csv", index=False)
ann = {}
with open(root + r"\all.genes.marker.stat.xls.xls", encoding="utf-8", errors="replace") as fh:
    fh.readline()
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 24: ann[p[1]] = p[14][:52]
print("\n=== 偏向【成熟效应血细胞(亚群1)】命运的基因 ===")
for _, r in res.tail(12).iloc[::-1].iterrows(): print("   %-20s rho=%+.2f | %s" % (r.gene, r.rho, ann.get(r.gene,"")))
print("\n=== 偏向【亚群0】命运的基因 ===")
for _, r in res.head(12).iterrows(): print("   %-20s rho=%+.2f | %s" % (r.gene, r.rho, ann.get(r.gene,"")))
A.write(out + r"\adata_hemocyte_lineage.h5ad")
A.obs[["sub","group","phase","S_score","G2M_score","dpt_pseudotime","fate_to_0","fate_to_1","fate_bias"]].to_csv(out + r"\fate_probabilities.csv")
B = A[A.obs_names[R]].copy()
B.obsm["X_umap"] = sc.read_h5ad(out + r"\adata_hemocyte_lineage.h5ad").obsm.get("X_umap", None) if False else None
sc.settings.figdir = fig
import scanpy
sc.pl.umap(A, color=["sub","phase","dpt_pseudotime","fate_to_1","fate_to_0"], show=False, save="_lineage.png", ncols=3)
print("\nfigures:", sorted(os.listdir(fig)))
