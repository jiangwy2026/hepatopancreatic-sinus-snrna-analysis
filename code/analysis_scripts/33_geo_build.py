import sys, io, os, glob, re
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 250); sc.settings.verbosity = 0
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"; GEO = out + r"\geo_ref"
num = lambda s: re.sub(r"^(ncbi_|LOC|gene-)", "", str(s))
files = sorted(glob.glob(os.path.join(GEO, "GSM*_filtered_feature_bc_matrix.h5")))
parts = []
for f in files:
    name = os.path.basename(f).split("_")[1]
    a = sc.read_10x_h5(f)
    a.var_names = [num(x) for x in a.var_names]
    a.var_names_make_unique()
    a.obs["sample"] = name
    a.obs_names = [name + "_" + b for b in a.obs_names]
    parts.append(a)
R = ad.concat(parts, join="inner", index_unique=None)
print("合并参考:", R.shape)
mine = pd.read_csv(root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\1.Expression\expressions\blank\features.tsv.gz", sep="\t", header=None)[0].astype(str).map(num).tolist()
print("基因 ID 归一化后交集: %d (占你的数据 %.1f%% | 占参考 %.1f%%)" % (len(set(mine) & set(R.var_names)), len(set(mine) & set(R.var_names))/len(set(mine))*100, len(set(mine) & set(R.var_names))/len(set(R.var_names))*100))
R.layers["counts"] = R.X.copy()
sc.pp.calculate_qc_metrics(R, percent_top=None, log1p=False, inplace=True)
print("\n参考数据 QC: 每细胞基因数分位 25%%=%.0f 50%%=%.0f 75%%=%.0f" % tuple(np.percentile(R.obs.n_genes_by_counts, [25,50,75])))
print("每细胞 counts 分位 25%%=%.0f 50%%=%.0f 75%%=%.0f" % tuple(np.percentile(R.obs.total_counts, [25,50,75])))
print("n_genes<200 的细胞: %d (%.1f%%)" % ((R.obs.n_genes_by_counts<200).sum(), (R.obs.n_genes_by_counts<200).mean()*100))
R.write(out + r"\geo_GSE306674_raw.h5ad")
print("saved geo_GSE306674_raw.h5ad")
