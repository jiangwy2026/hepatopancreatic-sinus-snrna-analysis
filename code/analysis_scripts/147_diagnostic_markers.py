# -*- coding: utf-8 -*-
"""
147_diagnostic_markers.py
判断"血细胞群匹配到肝胰腺类型"是真信号还是巧合：
1) 每套参考的 marker 数量（证据量）
2) 竞争身份的经典 marker 在各群的检出率（内皮/E细胞/角质细胞 vs 血细胞）
"""
import io, sys, os
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"
CL = ["0","1","2","3","4","5","6","7","8","9"]
OUR = {"0":"MH1","1":"GH","2":"PH","3":"MH2","4":"Epi","5":"F","6":"R","7":"B","8":"NEU","9":"NEU"}
DIAG = {
 "endothelial_canonical": ["CDH5","PECAM1","KDR","VWF","EMCN","CLDN5","CD34","SPON1"],
 "endothelial_ref(2 genes)": ["Igfbp7","SPON1"],
 "PH_progenitor": ["ccnd2","Grik2","PPAF2","SSPO"],
 "TGM(transglutaminase)": ["LOC113823945"],
 "Ecell_ref(2 genes)": ["Dsp1","eIF-5A"],
 "ubiquitous_control": ["Act5C","eef1a"],
 "granulocyte_ref": ["PEN-3","LOC113801825"],
 "phagocyte_ref": ["PRDX5","Nlrp3","ALF"],
}
A = sc.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))
A.X = A.layers["counts"].copy(); sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
sym = A.var["sym"].astype(str).values; ids = A.var_names.astype(str).values
tg = np.asarray(A.layers["counts"].sum(0)).ravel()
best = {}
for i, s in enumerate(sym):
    if s not in best or tg[i] > tg[best[s]]: best[s] = i
cl = A.obs["vendor_cluster"].astype(str).values
rows = []
for label, gl in DIAG.items():
    for g in gl:
        if g not in best:
            rows.append({"set": label, "gene": g, "note": "not found"}); continue
        v = np.asarray(A[:, best[g]].X.todense()).ravel()
        row = {"set": label, "gene": g, "note": ""}
        for c in CL: row["%s(C%s)" % (OUR[c], c)] = round(100 * (v[cl == c] > 0).mean(), 1)
        row["all"] = round(100 * (v > 0).mean(), 1)
        rows.append(row)
T = pd.DataFrame(rows)
cols = ["set","gene","all"] + ["%s(C%s)" % (OUR[c], c) for c in CL]
print("=== 诊断 marker 检出率（%）===")
print(T[cols].fillna("").to_string(index=False))
T.to_csv(os.path.join(OUT, "diagnostic_marker_pct.csv"), index=False)
print("\n=== 参考 marker 数量对比 ===")
print("  eLife 血细胞 6 亚型: PH1 746 / PH2 680 / MH1 626 / MH2 499 / GH2 172 / GH1 120（本数据匹配率 91.7–98.3%）")
print("  肝胰腺 12 类: F 4 / R 3 / B 3 / E 2 / M 2 / endothelial 2 / granulocyte 2 / phagocyte 3（每类 2–4 个）")
