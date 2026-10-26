import anndata as ad, os, pandas as pd
root = r"D:\桌面\single cell\analysis_output"
for f in ["adata_C2_within.h5ad","adata_hemocyte_traj.h5ad","adata_velocity_processed.h5ad"]:
    A = ad.read_h5ad(os.path.join(root,f))
    print("\n==",f,A.shape)
    print("obs",list(A.obs.columns))
    print("obsm",list(A.obsm.keys()))
    print("layers",list(A.layers.keys()))
    print(A.obs.head(3).to_string())
    for c in ["pt_bin","dpt_pseudotime","group","sample","phase","S_score","G2M_score"]:
        if c in A.obs: print(c, A.obs[c].value_counts(dropna=False).head(8).to_dict() if A.obs[c].dtype.name in ("category","object") else A.obs[c].describe().to_dict())
