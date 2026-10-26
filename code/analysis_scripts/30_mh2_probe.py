import sys, io, os
import numpy as np, pandas as pd, scanpy as sc
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 280)
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
A = sc.read_h5ad(out + r"\adata_soupx.h5ad")
A.X = A.layers["counts"].copy(); sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
sym2i = {}
for i, s in zip(A.var_names, A.var["sym"]): sym2i.setdefault(str(s), i)
probe = ["Nlrp3","Casp1","NLRP3","CASP1","VEGF3","VEGF5","LYZ1","Lyz1","CHIT1","ProPO","PEN-3","TGM1","Fli1"]
print("=== 论文招牌 marker 是否存在 ===")
found = {}
for g in probe:
    if g in sym2i: found[g] = sym2i[g]; print("  %-8s 存在" % g)
    else: print("  %-8s 不存在" % g)
# 模糊搜索 NLRP/Casp
for key in ("NLRP","Casp","VEGF","Lyz","ProPO"):
    hits = [s for s in sym2i if str(s).lower().startswith(key.lower())][:8]
    print("  含 '%s' 的基因: %s" % (key, hits))
if found:
    rows = {}
    for g, i in found.items():
        v = np.asarray(A[:, i].X.todense()).ravel()
        d = pd.DataFrame({"clu": A.obs.vendor_cluster.astype(str).values, "x": v})
        rows[g] = d.groupby("clu")["x"].mean().round(2)
    T = pd.DataFrame(rows).T
    print("\n=== 各群平均表达 ===")
    print(T.to_string())
