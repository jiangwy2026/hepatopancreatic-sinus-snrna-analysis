import sys, io
import numpy as np, pandas as pd, scanpy as sc
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 280)
out = r"D:\桌面\single cell\analysis_output"
A = sc.read_h5ad(out + r"\adata_soupx.h5ad")
A.X = A.layers["counts"].copy(); sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
sym2i = {}
for i, s in zip(A.var_names, A.var["sym"]): sym2i.setdefault(str(s), i)
for g in ["casp1-a","Nlrp3","Lyz2","CHIT1"]:
    if g not in sym2i: print(g, "不存在"); continue
    i = sym2i[g]; v = np.asarray(A[:, i].X.todense()).ravel()
    d = pd.DataFrame({"clu": A.obs.vendor_cluster.astype(str).values, "x": v})
    m = d.groupby("clu")["x"].mean().round(3); det = d.assign(p=v>0).groupby("clu")["p"].mean().mul(100).round(0)
    print("\n%s | 平均表达: %s" % (g, m.to_dict()))
    print("%s | 检出率%%: %s" % (" "*len(g), det.to_dict()))
