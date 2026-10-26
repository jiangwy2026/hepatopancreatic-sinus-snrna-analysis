import io, sys, os
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
from scipy.stats import rankdata, mannwhitneyu
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"
A = sc.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))
A.X = A.layers["counts"].copy(); sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
cl = A.obs["vendor_cluster"].astype(str).values
for gid in ("ncbi_113814652", "ncbi_113823802", "ncbi_113823945"):
    j = list(A.var_names.astype(str)).index(gid)
    x = np.asarray(A[:, j].X.todense()).ravel()
    m = cl == "2"; npos = m.sum(); nneg = len(x) - npos
    r = rankdata(x)
    auc_formula = (r[m].sum() - npos*(npos+1)/2)/(npos*nneg)
    u = mannwhitneyu(x[m], x[~m], alternative="two-sided").statistic / (npos*nneg)
    print("%s | pct PH %.1f%% vs rest %.1f%% | AUC formula %.3f | AUC via MannWhitney %.3f" % (
        A.var["sym"].values[j], 100*(x[m]>0).mean(), 100*(x[~m]>0).mean(), auc_formula, u))
