import sys, io, os, time
import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, scanpy as sc
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sc.settings.verbosity = 0
try:
    import openTSNE; print("openTSNE 已安装:", openTSNE.__version__)
except Exception:
    print("openTSNE 未安装 -> 将用 sklearn 的 t-SNE")
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
A = sc.read_h5ad(out + r"\adata_soupx_hvg.h5ad")
print("对象:", A.shape, "| obsm:", list(A.obsm.keys()))
t0 = time.time()
sc.tl.tsne(A, use_rep="X_pca_harmony", n_pcs=30, perplexity=30, random_state=0, n_jobs=8)
print("t-SNE 完成 (%.0fs)" % (time.time()-t0))
print("坐标范围: x %.1f~%.1f | y %.1f~%.1f" % (A.obsm["X_tsne"][:,0].min(), A.obsm["X_tsne"][:,0].max(),
                                               A.obsm["X_tsne"][:,1].min(), A.obsm["X_tsne"][:,1].max()))
A.write(out + r"\adata_soupx_hvg_tsne.h5ad")
print("saved adata_soupx_hvg_tsne.h5ad")
