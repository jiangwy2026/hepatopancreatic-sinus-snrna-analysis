import anndata as ad, os, pandas as pd
O=r"D:\桌面\single cell\analysis_output"
for f in ['elife_adata.h5ad','geo_GSE306674_processed.h5ad','adata_soupx.h5ad']:
 A=ad.read_h5ad(os.path.join(O,f))
 print('\n==',f,A.shape)
 print('obs',list(A.obs.columns))
 print('var',list(A.var.columns))
 print('obsm',list(A.obsm.keys()))
 print('layers',list(A.layers.keys()))
 print(A.obs.head(3).to_string())
 for c in A.obs.columns:
  if A.obs[c].dtype.name in ('category','object') and A.obs[c].nunique()<30:
   print(c,A.obs[c].value_counts(dropna=False).head(20).to_dict())
