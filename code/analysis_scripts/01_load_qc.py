import gzip, os, sys, io
import numpy as np, pandas as pd, scipy.io, scipy.sparse as sp
import anndata as ad, scanpy as sc
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

root = r"D:\桌面\single cell"
base = root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result"
expr = base + r"\1.Expression\expressions"
out  = root + r"\analysis_output"
MT = ["ND1","ND2","ND3","ND4","ND4L","ND5","COX1","COX2","COX3","ATP6","ATP8","CYTB"]

feats = pd.read_csv(expr + r"\blank\features.tsv.gz", sep="\t", header=None, names=["gid","sym","type"])

def load_sample(s):
    M = scipy.io.mmread(expr + "\\" + s + r"\matrix.mtx.gz").tocsr()
    bars = pd.read_csv(expr + "\\" + s + r"\barcodes.tsv.gz", sep="\t", header=None, names=["bc"])
    A = ad.AnnData(X=sp.csr_matrix(M.T))
    A.var_names = feats["gid"].values
    A.obs_names = [s + "_" + b.rsplit("-",1)[0] for b in bars["bc"]]
    A.obs["sample"] = s
    print(f"[{s}] {A.shape[0]} cells x {A.shape[1]} genes, nnz={A.X.nnz:,}")
    return A

A = ad.concat([load_sample("blank"), load_sample("treated")], join="outer", index_unique=None)
A.var["symbol"] = feats["sym"].values
A.var["mt"] = A.var["symbol"].isin(MT).values
A.layers["counts"] = A.X.copy()
print("cells:", A.shape[0], "genes:", A.shape[1], "MT genes:", int(A.var["mt"].sum()))
sc.pp.calculate_qc_metrics(A, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
A.write(out + r"\adata_raw.h5ad")
print("\n=== 原始数据 QC 概览 ===")
print(A.obs[["total_counts","n_genes_by_counts","pct_counts_mt"]].describe(percentiles=[.25,.5,.75,.95]).round(2).to_string())
print("\n=== 分样本中位数 ===")
print(A.obs.groupby("sample")[["total_counts","n_genes_by_counts","pct_counts_mt"]].median().round(2).to_string())
print("\nsaved", out + r"\adata_raw.h5ad")
