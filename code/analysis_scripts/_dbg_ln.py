import io, sys, os
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
OUT = r"D:\桌面\single cell\analysis_output"
raw = ad.read_h5ad(os.path.join(OUT, "adata_raw.h5ad"))
X = raw.layers["counts"].tocsr()
print("X dtype:", X.dtype, "shape:", X.shape, "nnz:", X.nnz)
tot = np.asarray(X.sum(1)).ravel()
norm = sp.diags(1.0 / np.maximum(tot, 1)) @ X
print("norm dtype:", norm.dtype, "nnz:", norm.nnz, "data[:5]:", norm.data[:5])
ln = norm.copy(); ln.data = np.log1p(ln.data)
gp = {g: i for i, g in enumerate(raw.var_names.astype(str))}
j = gp["ncbi_113812814"]
print("PH marker 列 j=%d 前 10 个细胞的 raw/norm/ln 值:" % j)
print("  raw :", X[:10, j].toarray().ravel())
print("  norm:", norm[:10, j].toarray().ravel())
print("  ln  :", ln[:10, j].toarray().ravel())
sub = ln[:50][:, [gp["ncbi_113812814"], gp["ncbi_113802728"]]]
print("  ln[:50][:,2cols].mean(1)[:5]:", np.asarray(sub.mean(1)).ravel()[:5])
