import sys, io, os, time
import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, scanpy as sc
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sc.settings.verbosity = 0
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
H = sc.read_h5ad(out + r"\adata_hemocyte_traj.h5ad")
base = {"2":"PH","1":"GH","0":"MH1"}
H.obs["cell"] = H.obs.vendor_cluster.astype(str).map(base)
H = H[H.obs["cell"].notna()].copy()
H.obs["cell"] = pd.Categorical(H.obs["cell"], categories=["PH","GH","MH1"])
print("血细胞:", H.shape, dict(H.obs["cell"].value_counts()))
if "X_umap" not in H.obsm:
    t0 = time.time()
    sc.pp.neighbors(H, use_rep="X_pca_harmony", n_neighbors=15, n_pcs=30)
    sc.tl.umap(H, random_state=0)
    print("UMAP 完成 (%.0fs)" % (time.time()-t0))
H.write(out + r"\adata_hemocyte_fig4.h5ad")
print("saved adata_hemocyte_fig4.h5ad")
