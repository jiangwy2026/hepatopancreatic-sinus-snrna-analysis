import sys, io, os
import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, scanpy as sc
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 250); sc.settings.verbosity = 0
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
F = sc.read_h5ad(out + r"\adata_soupx.h5ad")
H2 = F[F.obs.vendor_cluster.astype(str) == "2"].copy()
H2.layers["counts"] = H2.layers["counts"].copy(); H2.X = H2.layers["counts"].copy()
sc.pp.normalize_total(H2, target_sum=1e4); sc.pp.log1p(H2)
mean = np.asarray(H2.X.mean(0)).ravel(); sq = np.asarray(H2.X.multiply(H2.X).mean(0)).ravel()
var = sq - mean**2
hv = np.argsort(-var)[:2000]
Hb = H2[:, hv].copy(); Hb.X = Hb.layers["counts"].copy()
sc.pp.normalize_total(Hb, target_sum=1e4); sc.pp.log1p(Hb); sc.pp.scale(Hb, max_value=10)
sc.tl.pca(Hb, n_comps=30, svd_solver="arpack")
sc.pp.neighbors(Hb, n_neighbors=15, n_pcs=20)
lab = "c2sub"
sc.tl.leiden(Hb, resolution=0.3, key_added=lab, flavor="igraph", n_iterations=2, directed=False)
SC = ["PCNA","MCM6","CDC6","ORC1","POLA2","PRIM1","FEN1","TYMS","CCNB3","AURKA","TOP2A","TPX2","KIF20A","KIF23","RACGAP1","ASPM","Pka-C1","Patronin"]
PROG = ["LOC113823945","ccnd2","Igfbp7","Grik2","SSP0","slo","mbl","Zasp52"]
sym2i = {}
for i, s in zip(Hb.var_names, Hb.var["sym"].astype(str)): sym2i.setdefault(s, i)
res = {}
for nm, gl in (("周期基因", SC), ("祖细胞marker", PROG)):
    ids = [sym2i[g] for g in gl if g in sym2i]
    res[nm] = pd.Series(np.asarray(Hb[:, ids].X.mean(1)).ravel(), index=Hb.obs_names).groupby(Hb.obs[lab].values).mean()
T = pd.DataFrame(res); T["n"] = Hb.obs[lab].value_counts().sort_index()
T["treated占比%"] = (Hb.obs.groupby(lab)["group"].apply(lambda s:(s=="treated").mean()*100)).round(1)
print("=== C2 内部亚结构（res 0.3, %d 个亚群）===" % Hb.obs[lab].nunique()); print(T.round(3).to_string())
sc.tl.rank_genes_groups(Hb, lab, method="wilcoxon", pts=True)
mk = pd.concat([sc.get.rank_genes_groups_df(Hb, group=g).assign(sub=g) for g in Hb.obs[lab].cat.categories], ignore_index=True)
mk["gene"] = mk["names"].map(dict(zip(Hb.var_names, Hb.var["sym"].astype(str)))).fillna(mk["names"]).astype(str)
sig = mk[(mk.pvals_adj<0.05)&(mk.logfoldchanges>0.5)]
print("\n=== 各亚群 top marker ===")
for g in Hb.obs[lab].cat.categories:
    d = sig[sig["sub"]==g].sort_values("scores",ascending=False).head(8)
    print("  %-3s (n=%3d): %s" % (g, (Hb.obs[lab]==g).sum(), ", ".join(d.gene.tolist())))
