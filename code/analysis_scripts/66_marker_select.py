import sys, io, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np, pandas as pd, scanpy as sc
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 250); sc.settings.verbosity = 0
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"; figdir = out + r"\figures_pub"
A = sc.read_h5ad(out + r"\adata_soupx.h5ad")
base = {"2":"PH","1":"GH","0":"MH1","3":"MH2","4":"Epi","5":"F","6":"R","7":"B","8":"NE","9":"Neu"}
A.obs["ct"] = A.obs.vendor_cluster.astype(str).map(base).replace({"NE":"NEU","Neu":"NEU"})
ORD = ["PH","GH","MH1","MH2","Epi","F","R","B","NEU"]
A.obs["ct"] = pd.Categorical(A.obs["ct"], categories=ORD)
A.X = A.layers["counts"].copy(); sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
mt = {"ND1","ND2","ND3","ND4","ND4L","ND5","COX1","COX2","COX3","ATP6","ATP8","CYTB"}
symv = A.var["sym"].astype(str).values
A.var["is_mt"] = np.isin(symv, list(mt))
sc.tl.rank_genes_groups(A, "ct", method="wilcoxon", pts=True)
mk = pd.concat([sc.get.rank_genes_groups_df(A, group=g).assign(grp=g) for g in ORD], ignore_index=True)
mk["gene"] = mk["names"].map(dict(zip(A.var_names, symv))).fillna(mk["names"]).astype(str)
mk = mk.merge(A.var[["is_mt"]], left_on="names", right_index=True, how="left")
mk["named"] = ~mk["gene"].str.startswith(("LOC","gene-LOC"))
mk.to_csv(out + r"\markers_9groups.csv", index=False)
# 每群取 5 个：排除 MT；优先 padj<0.05 & log2FC>0.5 且有正式符号；不足则放宽
picked, seen = [], set()
for g in ORD:
    d = mk[(mk.grp == g) & (~mk.is_mt.fillna(False))].copy()
    d1 = d[(d.pvals_adj < 0.05) & (d.logfoldchanges > 0.5)]
    for cand in (d1.sort_values(["named","scores"], ascending=[False,False]),
                 d[d.pvals_adj < 0.05].sort_values(["named","scores"], ascending=[False,False]),
                 d.sort_values(["named","scores"], ascending=[False,False])):
        for _, r in cand.iterrows():
            if len([x for x in picked if x["grp"] == g]) >= 5: break
            if r.gene in seen: continue
            picked.append({"grp": g, "gene": r.gene, "log2FC": r.logfoldchanges, "q": r.pvals_adj})
            seen.add(r.gene)
        if len([x for x in picked if x["grp"] == g]) >= 5: break
P = pd.DataFrame(picked)
print("选出 %d 个 marker" % len(P))
print(P.groupby("grp", sort=False).gene.apply(lambda s: ", ".join(s)).to_string())
# 均值矩阵（z-score across groups）
sym2i = {}
for i, s in zip(A.var_names, symv): sym2i.setdefault(s, i)
M = np.zeros((len(P), len(ORD)))
for r_i, g in enumerate(P.gene):
    j = sym2i[g]
    for c_i, gg in enumerate(ORD):
        m = (A.obs["ct"] == gg).values
        M[r_i, c_i] = float(np.asarray(A[m][:, j].X.mean()).ravel()[0])
Z = (M - M.mean(1, keepdims=True)) / (M.std(1, keepdims=True) + 1e-9)
pd.DataFrame(Z, index=P.gene, columns=ORD).to_csv(out + r"\marker_heatmap_z.csv")
np.save(out + r"\_tmp_heatmap_Z.npy", Z); P.to_csv(out + r"\marker_heatmap_selected.csv", index=False)
print("\nz-score 范围: %.2f ~ %.2f" % (Z.min(), Z.max()))
