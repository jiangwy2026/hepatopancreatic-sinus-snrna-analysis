import anndata as ad, os
root = r"D:\桌面\single cell\analysis_output"
files = ["adata_hemocyte_traj_geodesic.h5ad","adata_hemocyte_traj_harmony-nn30.h5ad","adata_velocity_processed.h5ad","adata_soupx_hvg_tsne.h5ad","adata_soupx.h5ad"]
for f in files:
    p = os.path.join(root, f)
    A = ad.read_h5ad(p)
    print("\n==", f, A.shape)
    print("obs:", list(A.obs.columns))
    print("obsm:", list(A.obsm.keys()))
    print("uns:", list(A.uns.keys()))
    print("layers:", list(A.layers.keys()))
    print("obsp:", list(A.obsp.keys()))
    print(A.obs.head(3).to_string())
