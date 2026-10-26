# -*- coding: utf-8 -*-
"""160_C0_metabolic_scores.py - 核糖体/氧化磷酸化/TCA 程序在各群的相对水平（验证 C0 是'高代谢生物合成过渡态'）"""
import io, sys, os, re
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"
OUR = {"0":"MH1(C0)","1":"GH","2":"PH","3":"MH2","4":"Epi","5":"F","6":"R","7":"B","8":"Neu8","9":"Neu9"}
A = sc.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))
A.X = A.layers["counts"].copy(); sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
sym = A.var["sym"].astype(str).values
cl = A.obs["vendor_cluster"].astype(str).values
def genes(pat): return [i for i, s in enumerate(sym) if re.match(pat, s, re.I)]
PROG = {
 "Ribosome (RpL/RpS)": genes(r"^(Rp[LS]|rpl-|rps-|RPS|RPL)"),
 "Translation factors (eef/eIF)": genes(r"^(eef|eIF)"),
 "OXPHOS complex I (NDUF/ND)": [i for i, s in enumerate(sym) if re.match(r"^(NDUF|ND\d)", s, re.I)],
 "Cytochrome b/c (CYTB/COX)": [i for i, s in enumerate(sym) if re.match(r"^(CYTB|COX\d)", s, re.I)],
 "TCA cycle": [i for i, s in enumerate(sym) if re.match(r"^(SDH|MDH|CS|ACO|IDH|FH|OGDH|DLST|SUCL)", s, re.I)],
 "Proteasome": genes(r"^(PSM)"),
}
rows = {}
for name, gi in PROG.items():
    if not gi: continue
    vals = [float(np.asarray(A[cl == c][:, gi].X.mean())) for c in OUR]
    rows[name] = dict(zip(OUR.values(), np.round(vals, 3)))
T = pd.DataFrame(rows).T
print("=== 各群的代谢/生物合成程序（平均 log 表达）===")
print(T.to_string())
T.to_csv(os.path.join(OUT, "C0_metabolic_programs.csv"))
print("\n各程序最高的群:")
for p in T.index: print("  %-32s -> %s (%.3f)" % (p, T.loc[p].idxmax(), T.loc[p].max()))
