import sys, io, os, gzip, pickle
import numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad, scanpy as sc
from sklearn.metrics import adjusted_rand_score as ari
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
FIL = root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\1.Expression\expressions"
feats = pd.read_csv(FIL + r"\blank\features.tsv.gz", sep="\t", header=None, names=["gid","sym","type"])
store = pickle.load(open(out + r"\_tmp_soupx_store.pkl","rb"))
parts = []
for s in ("blank","treated"):
    C = store[s]["C"]; cells = store[s]["cells"]; clu = store[s]["clu"]
    A = ad.AnnData(X=C)
    A.obs_names = cells; A.var_names = feats["gid"].astype(str).values
    A.obs["sample"] = s; A.obs["group"] = s; A.obs["vendor_cluster"] = pd.Categorical(clu.astype(str))
    parts.append(A)
A = ad.concat(parts, join="outer", index_unique=None)
A.var["sym"] = feats["sym"].astype(str).values
A.var["mt"] = A.var["sym"].isin(["ND1","ND2","ND3","ND4","ND4L","ND5","COX1","COX2","COX3","ATP6","ATP8","CYTB"]).values
A.layers["counts"] = A.X.copy()
sc.pp.calculate_qc_metrics(A, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
print("校正后对象:", A.shape)
print("每细胞 counts 中位 %.0f | 基因数中位 %.0f | pct_mt 中位 %.2f" % (np.median(A.obs.total_counts), np.median(A.obs.n_genes_by_counts), np.median(A.obs.pct_counts_mt)))
A.write(out + r"\adata_soupx.h5ad")
print("saved adata_soupx.h5ad")
