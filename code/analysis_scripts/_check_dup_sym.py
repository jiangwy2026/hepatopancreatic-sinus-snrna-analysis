import os, numpy as np, pandas as pd, scanpy as sc
OUT = r"D:\桌面\single cell\analysis_output"
A = sc.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))
sym = A.var["sym"].astype(str).values
vc = pd.Series(sym).value_counts()
dup = vc[vc > 1]
print("总基因数:", len(sym), "| 唯一符号数:", vc.size, "| 重复符号:", len(dup), "| 涉及基因数:", int(dup.sum()))
print("重复最多的符号:", dup.head(8).to_dict())
for g in ("GRIK2", "Pxt", "ccnd2", "PEN-3", "A2M"):
    idx = np.where(sym == g)[0]
    tot = np.asarray(A.layers["counts"][:, idx].sum(0)).ravel()
    print("%-6s 出现 %d 次 | 各基因ID/总counts: %s" % (g, len(idx), [(A.var_names[i], int(t)) for i, t in zip(idx, tot)]))
mk = pd.read_csv(os.path.join(OUT, "blood4_dotplot_genes.csv"))
print()
print("marker 表里记录的 gid 与按符号取到的 gid 是否一致：")
sym2first = {}
for i, s in enumerate(sym): sym2first.setdefault(s, A.var_names[i])
bad = [(r.sym, r.gid, sym2first[r.sym]) for r in mk.itertuples() if str(r.gid) != str(sym2first[r.sym])]
print("不一致的数量:", len(bad), "of", len(mk))
for b in bad[:12]: print("   ", b)
