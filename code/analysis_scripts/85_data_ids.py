# -*- coding: utf-8 -*-
import gzip, re, collections
base = r"D:\桌面\single cell\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\1.Expression\expressions"
for s in ("blank", "treated"):
    with gzip.open(base + "\\" + s + "\\features.tsv.gz", "rt") as f:
        rows = [l.rstrip("\n").split("\t") for l in f]
    with gzip.open(base + "\\" + s + "\\barcodes.tsv.gz", "rt") as f:
        bc = [l.strip() for l in f]
    print(s, "n_features:", len(rows), "cols:", len(rows[0]))
    print("   first 3 features:", rows[:3])
    print("   n_barcodes:", len(bc), "first 3:", bc[:3])
    ids = [r[0] for r in rows]; names = [r[1] if len(r) > 1 else "" for r in rows]
    print("   any 'ncbi_' in col0:", sum(x.startswith("ncbi_") for x in ids), "| sample col1:", names[:3])
print()
gtf = r"E:\workspace\genome2\genome\GENOME_GTF\GCF_003789085.1.gtf.gz"
n = 0
with gzip.open(gtf, "rt") as f:
    for line in f:
        if line.startswith("#"): continue
        p = line.split("\t")
        if len(p) > 8 and p[2] == "gene":
            print("GTF gene line:", p[8].strip()[:160]); n += 1
            if n >= 3: break
print()
import anndata as ad
a = ad.read_h5ad(r"D:\桌面\single cell\analysis_output\adata_soupx.h5ad", backed="r")
print("adata:", a.shape, "obs_names[:3]:", list(a.obs_names[:3]), "var_names[:3]:", list(a.var_names[:3]))
print("obs columns:", list(a.obs.columns))
if "sample" in a.obs.columns:
    print("sample counts:", a.obs["sample"].value_counts().to_dict())
