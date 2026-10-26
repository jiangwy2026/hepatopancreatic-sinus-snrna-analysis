import sys, io, os, re
import numpy as np, pandas as pd, scanpy as sc
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 300)
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
num = lambda s: re.sub(r"^(ncbi_|LOC|gene-)", "", str(s))
mk = pd.read_csv(out + r"\geo_ref_cluster_markers_sym.csv")
A = sc.read_h5ad(out + r"\adata_soupx.h5ad")
A.X = A.layers["counts"].copy(); sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
A.var_names = [num(x) for x in A.var_names]; A.var_names_make_unique()
myg = set(A.var_names)
mine_clu = A.obs.vendor_cluster.astype(str).values
# 参考各群 top50 marker（限定共享基因）
sets = {}
for g in sorted(mk["cl"].unique(), key=lambda x: int(x)):
    d = mk[mk["cl"] == g].sort_values("scores", ascending=False)
    ids = [num(x) for x in d["names"].astype(str).tolist()]
    ids = [i for i in dict.fromkeys(ids) if i in myg][:50]
    if len(ids) >= 10: sets["R%s" % g] = ids
print("可用于打分的参考群: %d 个" % len(sets))
res = {}
for lab, ids in sets.items():
    sub = A[:, ids]
    vals = {}
    for c in sorted(set(mine_clu), key=int):
        m = mine_clu == c
        vals[c] = float(np.asarray(sub.X[m].mean(0)).ravel().mean())
    res[lab] = pd.Series(vals)
S = pd.DataFrame(res).T
Z = S.sub(S.mean(axis=1), axis=0).div(S.std(axis=1) + 1e-9, axis=0)
print("\n=== 参考群 marker 在你各群的平均表达（行=参考群）===")
print(S.round(3).to_string())
print("\n=== 行内 z-score（每行最高 = 你数据中对应的群）===")
print(Z.round(2).to_string())
print("\n=== 你的每个群 匹配到哪个参考群 ===")
for c in Z.columns:
    o = Z[c].sort_values(ascending=False)
    print("  C%-3s -> %-4s (z=%+.2f)  次优 %-4s (z=%+.2f)" % (c, o.index[0], o.iloc[0], o.index[1], o.iloc[1]))
S.to_csv(out + r"\geo_label_transfer_scores.csv"); Z.to_csv(out + r"\geo_label_transfer_z.csv")
print("\nsaved geo_label_transfer_*.csv")
