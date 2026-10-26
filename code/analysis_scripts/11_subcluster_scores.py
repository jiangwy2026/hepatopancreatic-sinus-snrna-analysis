import sys, io, os
import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, scanpy as sc
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 250)
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
Hh = sc.read_h5ad(out + r"\adata_hemocyte_harmony.h5ad")
sub = Hh.obs["sub_0.3"].astype(str)
A = sc.read_h5ad(out + r"\adata_raw.h5ad"); A = A[A.obs_names.isin(set(Hh.obs_names))].copy()
A.X = A.layers["counts"].copy(); sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
A.obs["sub"] = sub.loc[A.obs_names].values
sym2id = {}
for i, s in zip(A.var_names, A.var["symbol"]): sym2id.setdefault(str(s), i)
sets = {
 "S_effector": ["A2M","PPAF1","PPAF2","Pxt","FCN2","LOC113803873","pigs","Ces2","Inx3"],
 "S_prolif":   ["ccnd2","LOC113823945","SSPO","slo","Grik2","Pka-C1","Patronin"],
 "S_amp":      ["PEN-3","LOC113801825","LOC113800363","LOC113813598","LOC113815775"],
 "S_hepato":   ["Prss1","Prss2","Try3","Cpa2","ACE","Oacyl","Mal-A2","LIPF"],
 "S_housekeep":["Act5C","eef1a","RpLP1","RpL10","LOC113830566","CYPA","LOC113822661"],
 "S_apod":     ["Apod","LOC113824714","LOC113728","LOC113824715"],
}
for k, gl in sets.items():
    ids = [sym2id[g] for g in gl if g in sym2id]
    miss = [g for g in gl if g not in sym2id]
    if miss: print("  %s 缺: %s" % (k, miss))
    sc.tl.score_genes(A, ids, score_name=k, ctrl_size=50)
cols = [k for k in sets]
q = pd.read_csv(out + r"\qc_per_cell_raw.csv").set_index("Cells")
A.obs["nCount"] = q.loc[A.obs_names,"nCount"].values
A.obs["nFeature"] = q.loc[A.obs_names,"nFeature"].values
S = A.obs.groupby("sub")[cols+["nCount","nFeature","pct_counts_mt"]].mean().round(3)
S["n_cells"] = A.obs.groupby("sub").size()
print("\n=== 亚群打分（平均）==="); print(S.to_string())
S.to_csv(out + r"\subcluster_scores.csv")
sig = pd.read_csv(out + r"\subcluster_markers_final.csv")
sig = sig[(sig.pvals_adj<0.05)&(sig.logfoldchanges>0.5)]
genes = []
for g in sorted(A.obs["sub"].unique()):
    genes += sig[sig["sub"]==g].sort_values("scores",ascending=False).head(8).gene.astype(str).tolist()
genes = [g for g in dict.fromkeys(genes) if g in sym2id]
A2 = A.copy(); A2.var_names = [str(x) for x in A2.var["symbol"].values]; A2.var_names_make_unique()
sc.settings.figdir = out + r"\figures_sub"
sc.pl.dotplot(A2, genes, groupby="sub", show=False, save="_markers_dotplot.png")
sc.pl.matrixplot(A2, genes, groupby="sub", show=False, save="_markers_matrix.png", cmap="Reds")
print("\nfigures:", sorted(os.listdir(out + r"\figures_sub")))
