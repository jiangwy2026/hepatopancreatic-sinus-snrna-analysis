import sys, io, os
import numpy as np, pandas as pd, scanpy as sc
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 280)
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"; REF = out + r"\elife_ref"
print("=== supp1 / supp2 的 Target_Cluster 取值与 marker 数 ===")
for f in ("elife-80127-supp1-v2.csv","elife-80127-supp2-v2.csv"):
    d = pd.read_csv(os.path.join(REF, f), encoding="utf-8-sig")
    print("\n%s  行数=%d" % (f, len(d)))
    print(d.groupby("Target_Cluster").size().to_string())
    print("  代表基因:", ", ".join(d["Gene ID"].astype(str).head(6).tolist()))
# 用 supp1/supp2 做群级打分
A = sc.read_h5ad(out + r"\adata_soupx.h5ad")
A.X = A.layers["counts"].copy(); sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
myg = set(A.var_names)
for f in ("elife-80127-supp1-v2.csv","elife-80127-supp2-v2.csv"):
    d = pd.read_csv(os.path.join(REF, f), encoding="utf-8-sig")
    d = d[d.Qvalue < 0.05]
    print("\n########## %s" % f)
    tab = {}
    for c, g in d.groupby("Target_Cluster"):
        ids = [x for x in g["Gene ID"].astype(str).unique() if x in myg]
        if len(ids) < 5: continue
        sub = A[:, ids]
        vals = {}
        for cl in sorted(A.obs.vendor_cluster.astype(str).unique(), key=int):
            m = (A.obs.vendor_cluster.astype(str) == cl).values
            vals[cl] = float(np.asarray(sub.X[m].mean(0)).ravel().mean())
        tab["cluster_%s(n=%d)" % (c, len(ids))] = pd.Series(vals)
    T = pd.DataFrame(tab)
    Z = T.sub(T.mean(0), axis=1).div(T.std(0) + 1e-9, axis=1)
    print("行内 z-score（每列最高=最匹配）:"); print(Z.round(2).to_string())
