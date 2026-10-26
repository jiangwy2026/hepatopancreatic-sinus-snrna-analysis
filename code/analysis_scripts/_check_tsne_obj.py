import anndata as ad, numpy as np
A = ad.read_h5ad(r"D:\桌面\single cell\analysis_output\adata_soupx_hvg_tsne.h5ad", backed="r")
print("tsne obj:", A.shape, "| obsm:", list(A.obsm.keys()), "| obs:", [c for c in A.obs.columns][:14])
G = ad.read_h5ad(r"D:\桌面\single cell\analysis_output\adata_hemocyte_traj_geodesic.h5ad", backed="r")
print("geodesic obj:", G.shape, "| has X_tsne:", "X_tsne" in A.obsm.keys())
names = set(map(str, A.obs_names))
print("geodesic cells present in tsne obj:", sum(1 for n in map(str, G.obs_names) if n in names), "of", G.n_obs)
if "phase" in A.obs.columns:
    print("phase counts:", A.obs["phase"].value_counts().to_dict())
