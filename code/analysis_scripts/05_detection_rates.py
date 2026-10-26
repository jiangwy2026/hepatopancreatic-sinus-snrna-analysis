import sys, io
import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, scanpy as sc
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
A = sc.read_h5ad(out + r"\adata_raw.h5ad")
cl = pd.read_csv(root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\2.Cluster\2.cluster\Cells.cluster.list.xls", sep="\t")
cl["Cells"] = cl["Cells"].astype(str)
A = A[A.obs_names.isin(set(cl["Cells"]))].copy()
A.obs["clu"] = cl.set_index("Cells").loc[A.obs_names, "Cluster"].astype(str).values
sym2id = {}
for i, s in zip(A.var_names, A.var["symbol"]): sym2id.setdefault(s, i)
genes = ["PEN-3","LOC113801825","LOC113800363","LOC113813598","Prss1","Try3","ACE","Oacyl","Mal-A2",
         "LOC113809538","A2M","LOC113823945","RASGRF1","LOC113803948","ccnd2"]
rows = {}
for g in genes:
    v = A[:, sym2id[g]].layers["counts"].toarray().ravel()
    d = pd.DataFrame({"clu": A.obs["clu"].values, "c": v})
    rows[g+"_pct"] = (d.groupby("clu")["c"].apply(lambda x: (x > 0).mean()*100)).round(0).astype(int)
    rows[g+"_mean"] = d[d.c > 0].groupby("clu")["c"].mean().round(1)
T = pd.DataFrame(rows).T
T.index.name = "gene"
pd.set_option("display.width", 250)
print("=== 检出率(%) 与阳性细胞平均 counts ===")
print(T.to_string())
print("\n=== 组间差异：每群 blank vs treated 占比变化 ===")
cm = pd.read_csv(out + r"\cellmeta_raw.csv")
t = pd.crosstab(cm.Cluster, cm.Groups)
t["pct_blank"] = (t["blank"]/t["blank"].sum()*100).round(1)
t["pct_treated"] = (t["treated"]/t["treated"].sum()*100).round(1)
t["fold_change"] = (t["pct_treated"]/t["pct_blank"]).round(2)
print(t.to_string())
