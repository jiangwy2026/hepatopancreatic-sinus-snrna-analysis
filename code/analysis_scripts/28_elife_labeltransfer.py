import sys, io, os
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 280)
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"; REF = out + r"\elife_ref"
maps = {"elife-80127-supp3-v2.csv":"PH1(原血细胞1)","elife-80127-supp4-v2.csv":"GH2(颗粒2)",
        "elife-80127-supp5-v2.csv":"MH2(单核2·巨噬样)","elife-80127-supp6-v2.csv":"PH2(原血细胞2)",
        "elife-80127-supp7-v2.csv":"GH1(颗粒1)","elife-80127-supp8-v2.csv":"MH1(单核1)"}
sets = {}
for f, lab in maps.items():
    d = pd.read_csv(os.path.join(REF, f), encoding="utf-8-sig")
    d = d[(d.Qvalue < 0.05)]
    ids = d["Gene ID"].astype(str).str.strip().unique().tolist()
    sets[lab] = ids
    print("%-22s %d 个 marker" % (lab, len(ids)))

A = sc.read_h5ad(out + r"\adata_soupx.h5ad")
A.X = A.layers["counts"].copy(); sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
mygenes = set(A.var_names)
print("\n=== marker 基因在你数据中的匹配率 ===")
for lab, ids in sets.items():
    hit = [g for g in ids if g in mygenes]
    print("  %-22s 匹配 %4d / %4d (%.1f%%)" % (lab, len(hit), len(ids), len(hit)/len(ids)*100))
# 计算每群平均表达
rows = {}
for lab, ids in sets.items():
    hit = [g for g in ids if g in mygenes]
    sub = A[:, hit]
    m = pd.Series(np.asarray(sub.X.mean(0)).ravel(), index=hit)
    per = A.obs.assign(_x=0)
    # 按群求均值
    expr = pd.DataFrame({"g": sub.var_names, "mean": np.asarray(sub.X.mean(0)).ravel()})
    # 逐群
    vals = {}
    for c in sorted(A.obs.vendor_cluster.astype(str).unique(), key=int):
        idx = (A.obs.vendor_cluster.astype(str) == c).values
        vals[c] = float(np.asarray(sub.X[idx].mean(0)).ravel().mean())
    rows[lab] = pd.Series(vals)
S = pd.DataFrame(rows).T
print("\n=== 交叉打分：eLife 亚型 marker 在你各群的【平均表达】 ===")
print(S.round(3).to_string())
Z = S.sub(S.mean(1), axis=0).div(S.std(1) + 1e-9, axis=0)
print("\n=== 同一矩阵的行内 z-score（每行最高=最匹配的群）===")
print(Z.round(2).to_string())
print("\n=== 每个 eLife 亚型在你的数据中匹配最好的群 ===")
for lab in Z.index:
    o = Z.loc[lab].sort_values(ascending=False)
    print("  %-22s -> C%-3s (z=%+.2f)  次优 C%s (z=%+.2f)" % (lab, o.index[0], o.iloc[0], o.index[1], o.iloc[1]))
print("\n=== 你的每个群最匹配哪个 eLife 亚型 ===")
for c in Z.columns:
    o = Z[c].sort_values(ascending=False)
    print("  C%-3s -> %-22s (z=%+.2f)  次优 %s (z=%+.2f)" % (c, o.index[0], o.iloc[0], o.index[1], o.iloc[1]))
S.to_csv(out + r"\elife_label_transfer_scores.csv")
Z.to_csv(out + r"\elife_label_transfer_z.csv")
print("\nsaved elife_label_transfer_scores.csv / _z.csv")
