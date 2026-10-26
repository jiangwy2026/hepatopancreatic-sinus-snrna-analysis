import anndata as ad, numpy as np, os, sys
from sklearn.manifold import TSNE
sys.stdout.reconfigure(encoding="utf-8")
OUT=r"D:\桌面\single cell\analysis_output"; STEP2=os.path.join(OUT,"step2"); STEP3=os.path.join(OUT,"step3")
a=ad.read_h5ad(os.path.join(STEP2,"step2_integrated.h5ad"))
Z=a.obsm["X_pca_harmony"]
print("running t-SNE", Z.shape, flush=True)
ts=TSNE(n_components=2, perplexity=30, init="pca", learning_rate="auto", max_iter=1000,
        random_state=0, n_jobs=-1).fit_transform(Z)
a.obsm["X_tsne"]=ts
a.write_h5ad(os.path.join(STEP2,"step2_integrated.h5ad"), compression="gzip")
b=ad.read_h5ad(os.path.join(STEP3,"step3_lyo_scores.h5ad"))
b.obsm["X_tsne"]=ts
b.write_h5ad(os.path.join(STEP3,"step3_lyo_scores.h5ad"), compression="gzip")
print("done", ts.shape)
