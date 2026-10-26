import anndata as ad, numpy as np, pandas as pd, scanpy as sc, os
root=r"D:\桌面\single cell\analysis_output"
A=ad.read_h5ad(os.path.join(root,"adata_C2_within.h5ad"))
print(A.obs.groupby(["pt_bin","group"]).size())
print(A.obs.groupby("pt_bin")["dpt_pseudotime"].agg(["min","max","median"]))
sym2gid={str(s):str(g) for s,g in zip(A.var["sym"].astype(str),A.var_names.astype(str))}
for g in ["ccnd2","LOC113823945","TGM1","Igfbp7"]:
 print(g, sym2gid.get(g))
