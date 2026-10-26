import sys, io, os
import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, scanpy as sc
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
Hh = sc.read_h5ad(out + r"\adata_hemocyte_harmony.h5ad")
A = sc.read_h5ad(out + r"\adata_raw.h5ad"); A = A[A.obs_names.isin(set(Hh.obs_names))].copy()
A.X = A.layers["counts"].copy(); sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
A.obs["sub"] = Hh.obs["sub_0.3"].astype(str).loc[A.obs_names].values
sig = pd.read_csv(out + r"\subcluster_markers_final.csv"); sig = sig[(sig.pvals_adj<0.05)&(sig.logfoldchanges>0.5)]
panel = ["A2M","PPAF2","Pxt","FCN2","PPAF1","pigs",
         "ccnd2","LOC113823945","SSPO","Grik2","slo","Igfbp7",
         "THBS4","Act5C","LOC113819543","LOC113830563","eef1a","RpLP1",
         "Apod","LOC113824714","LOC113824715","Crabp2","LOC113827839","LOC113828299"]
A2 = A.copy(); A2.var_names = [str(x) for x in A2.var["symbol"].values]; A2.var_names_make_unique()
panel = [g for g in panel if g in set(A2.var_names)]
sc.settings.figdir = out + r"\figures_sub"
sc.pl.dotplot(A2, panel, groupby="sub", show=False, save="_marker_panel_dotplot.png", figsize=(13,4), standard_scale="var")
sc.pl.matrixplot(A2, panel, groupby="sub", show=False, save="_marker_panel_matrix.png", figsize=(13,4), cmap="Reds", standard_scale="var")
sc.pl.stacked_violin(A2, panel, groupby="sub", show=False, save="_marker_panel_violin.png", figsize=(13,4), standard_scale="var")
print("figures:", sorted(os.listdir(out + r"\figures_sub")))
