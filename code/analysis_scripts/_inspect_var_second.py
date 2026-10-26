import anndata as ad, os, pandas as pd
O=r"D:\桌面\single cell\analysis_output"
for f in ['elife_adata.h5ad','geo_GSE306674_processed.h5ad']:
 A=ad.read_h5ad(os.path.join(O,f))
 print('\n',f)
 print(A.var.head(10).to_string())
 print('var_names',A.var_names[:10].tolist())
