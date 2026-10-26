import anndata as ad, numpy as np
A = ad.read_h5ad(r"D:\桌面\single cell\analysis_output\adata_soupx_hvg.h5ad", backed="r")
B = ad.read_h5ad(r"D:\桌面\single cell\analysis_output\adata_velocity_raw.h5ad", backed="r")
print("hvg:", A.shape, "| obs:", list(A.obs.columns), "| obsm:", list(A.obsm.keys()), "| layers:", list(A.layers.keys()))
print("velocity_raw:", B.shape, "| obs:", list(B.obs.columns), "| layers:", list(B.layers.keys()))
print("same cells:", bool((A.obs_names == B.obs_names).all()))
if "leiden" in A.obs.columns:
    print("leiden categories:", A.obs["leiden"].astype(str).value_counts().to_dict())
