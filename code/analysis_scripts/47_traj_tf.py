import sys, io, os
import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
from scipy.stats import rankdata
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 250); sc.settings.verbosity = 0
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
H = sc.read_h5ad(out + r"\adata_hemocyte_soupx.h5ad")
H.X = H.layers["counts"].copy(); sc.pp.normalize_total(H, target_sum=1e4); sc.pp.log1p(H)
# TF 注释
tf = {}
with open(root + r"\all.genes.marker.stat.xls.xls", encoding="utf-8", errors="replace") as fh:
    hd = fh.readline().rstrip("\n").split("\t")
    i_sym, i_tf = hd.index("Symbol"), hd.index("TF_family")
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) > i_tf and p[i_tf] not in ("-", ""):
            tf[p[i_sym]] = p[i_tf]
print("有 TF 家族注释的基因: %d" % len(tf))
from collections import Counter
print("  主要家族:", dict(Counter(tf.values()).most_common(8)))
# 拟时序相关
X = H.X; pt = H.obs["dpt_pseudotime"].values.astype(float)
rng = np.random.default_rng(0); si = rng.choice(len(pt), size=min(5000,len(pt)), replace=False)
ry = rankdata(pt[si]); ry = (ry-ry.mean())/(ry.std()+1e-12)
hv = np.argsort(-np.asarray(H.X[si].multiply(H.X[si]).mean(0)).ravel() - np.asarray(H.X[si].mean(0)).ravel()**2)[:5000]
Xs = X[si][:, hv].toarray().astype(np.float32)
rx = rankdata(Xs, axis=0); rx = (rx-rx.mean(0))/(rx.std(0)+1e-12)
rho = (rx*ry[:,None]).mean(0)
res = pd.DataFrame({"sym": H.var["sym"].values[hv].astype(str), "rho": rho})
res["TF"] = res["sym"].map(tf)
res = res.dropna(subset=["rho"]).sort_values("rho")
res.to_csv(out + r"\hemocyte_genes_vs_pseudotime.csv", index=False)
print("\n=== 拟时序【下降】= 祖细胞端特征 top20 ===")
for _,r in res.head(20).iterrows():
    print("   %-22s rho=%+.2f %s" % (r.sym, r.rho, ("[TF:%s]" % r.TF) if pd.notna(r.TF) else ""))
print("\n=== 拟时序【上升】= 成熟端特征 top20 ===")
for _,r in res.tail(20).iloc[::-1].iterrows():
    print("   %-22s rho=%+.2f %s" % (r.sym, r.rho, ("[TF:%s]" % r.TF) if pd.notna(r.TF) else ""))
# 只看 TF
tfres = res.dropna(subset=["TF"]).copy()
print("\n=== 转录因子随拟时序变化（全部 %d 个）— 下降端 top12（祖细胞驱动 TF）===" % len(tfres))
for _,r in tfres.head(12).iterrows(): print("   %-20s rho=%+.2f [%s]" % (r.sym, r.rho, r.TF))
print("=== 上升端 top12（成熟驱动 TF）===")
for _,r in tfres.tail(12).iloc[::-1].iterrows(): print("   %-20s rho=%+.2f [%s]" % (r.sym, r.rho, r.TF))
tfres.to_csv(out + r"\hemocyte_TFs_vs_pseudotime.csv", index=False)
# 关键标记基因随拟时序
panel = ["LOC113823945","ccnd2","Igfbp7","Grik2","SSPO","PEN-3","LOC113801825","A2M","Pxt","LOC113803873"]
print("\n=== 关键基因与拟时序的相关性 ===")
for g in panel:
    r = res[res.sym==g]
    if len(r): print("   %-20s rho=%+.2f" % (g, r.rho.iloc[0]))
    else: print("   %-20s (未进入 top5000 高变基因)" % g)
H.write(out + r"\adata_hemocyte_traj.h5ad")
print("\nsaved")
