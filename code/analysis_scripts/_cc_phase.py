# -*- coding: utf-8 -*-
import os, numpy as np, pandas as pd, scanpy as sc
OUT = r"D:\桌面\single cell\analysis_output"
ROWS = ["PH", "GH", "MH1", "MH2"]
S = ["PCNA","MCM6","CDC6","ORC1","ORC4","ORC5","POLA2","PRIM1","FEN1","TYMS","DHFR","HIRA","CHAF1A"]
G2M = ["CCNB3","AURKA","TOP2A","TPX2","ASPM","Pka-C1","KIF20A","KIF23","RACGAP1","NDC80","INCENP","PRC1","Patronin"]
A = sc.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))
base = {"2": "PH", "1": "GH", "0": "MH1", "3": "MH2"}
ct = A.obs["vendor_cluster"].astype(str).map(base)
H = A[ct.isin(ROWS).values].copy()
H.obs["ct"] = pd.Categorical(ct[ct.isin(ROWS)].values, categories=ROWS)
H.X = H.layers["counts"].copy(); sc.pp.normalize_total(H, target_sum=1e4); sc.pp.log1p(H)
sym2id = {}
for i, s in enumerate(H.var["sym"].astype(str).values):
    sym2id.setdefault(s, H.var_names[i])
s_ids = [sym2id[g] for g in S if g in sym2id]; g2m_ids = [sym2id[g] for g in G2M if g in sym2id]
print("S ids:", len(s_ids), "| G2M ids:", len(g2m_ids))
sc.tl.score_genes_cell_cycle(H, s_genes=s_ids, g2m_genes=g2m_ids)
comp = pd.crosstab(H.obs["ct"], H.obs["phase"], normalize="index").mul(100).round(1)
print("\nphase composition per group (%):")
print(comp.to_string())
print("\nmean scores:")
print(H.obs.groupby("ct", observed=True)[["S_score", "G2M_score"]].mean().round(4).to_string())
print("\ncell cycle gene detection (%% of cells with >0):")
gi = [H.var_names.get_loc(x) for x in s_ids + g2m_ids]
X = H[:, gi].X
X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
det = pd.DataFrame({"ct": H.obs["ct"].astype(str).values,
                    "ngene": (X > 0).sum(1)}).groupby("ct")["ngene"].agg(["mean", "median"])
print(det.round(2).to_string())
comp.to_csv(os.path.join(OUT, "cellcycle_phase_composition.csv"))
