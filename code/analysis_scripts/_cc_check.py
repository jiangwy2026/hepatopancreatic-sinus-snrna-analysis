# -*- coding: utf-8 -*-
import os, numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
OUT = r"D:\桌面\single cell\analysis_output"
ROWS = ["PH", "GH", "MH1", "MH2"]
CC = {"G1/S": ["PCNA","MCM6","CDC6","ORC1","ORC4","ORC5","POLA2","PRIM1","FEN1","TYMS","DHFR","HIRA","CHAF1A"],
      "G2":   ["CCNB3","AURKA","TOP2A","TPX2","ASPM","Pka-C1"],
      "M":    ["KIF20A","KIF23","RACGAP1","NDC80","INCENP","PRC1","Patronin"]}
A = sc.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))
base = {"2": "PH", "1": "GH", "0": "MH1", "3": "MH2"}
ct = A.obs["vendor_cluster"].astype(str).map(base)
H = A[ct.isin(ROWS).values].copy()
H.obs["ct"] = pd.Categorical(ct[ct.isin(ROWS)].values, categories=ROWS)
H.X = H.layers["counts"].copy(); sc.pp.normalize_total(H, target_sum=1e4); sc.pp.log1p(H)
sym2i = {}
for i, s in enumerate(H.var["sym"].astype(str).values): sym2i.setdefault(s, i)
for ph, gl in CC.items():
    genes = [g for g in gl if g in sym2i]
    X = H[:, [sym2i[g] for g in genes]].X
    X = np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X)
    pct = np.vstack([100*(X[(H.obs["ct"]==g).values] > 0).mean(0) for g in ROWS])
    mean = np.vstack([X[(H.obs["ct"]==g).values].mean(0) for g in ROWS])
    print("===", ph, "genes:", len(genes))
    print(pd.DataFrame(pct, index=ROWS, columns=genes).round(1).to_string())
    print("   mean %expr per row:", dict(zip(ROWS, pct.mean(1).round(1))))
    print("   mean log-expr per row:", dict(zip(ROWS, mean.mean(1).round(2))))
# independent: phase composition
sc.tl.score_genes_cell_cycle(H, s_genes=[g for g in CC["G1/S"] if g in sym2i],
                             g2m_genes=[g for g in CC["G2"]+CC["M"] if g in sym2i])
print("\nphase composition per group (%):")
comp = pd.crosstab(H.obs["ct"], H.obs["phase"], normalize="index").mul(100).round(1)
print(comp.to_string())
print("\nmean S_score / G2M_score:")
print(H.obs.groupby("ct", observed=True)[["S_score","G2M_score"]].mean().round(3).to_string())
