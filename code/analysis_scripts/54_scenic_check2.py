import sys, io, os, pandas as pd, numpy as np
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"; DB = out + r"\scenic_db"
import pyarrow.feather as feather
tbl = feather.read_table(os.path.join(DB, "dm6_v10_clust.genes_vs_motifs.rankings.feather"), memory_map=True)
print("=== gene_based 排序库 ===")
print("  %d 基因 x %d motif" % (tbl.num_rows, tbl.num_columns))
print("  基因名前 8:", [str(x) for x in tbl.column(0)[:8]])
print("  motif 前 3:", [str(x) for x in tbl.column_names[:3]])
feat = set(str(x) for x in tbl.column(0))
feats = pd.read_csv(root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\1.Expression\expressions\blank\features.tsv.gz", sep="\t", header=None, names=["gid","sym","type"])
mysym = set(feats["sym"].astype(str))
hit = mysym & feat
print("\n  我的基因符号 ∩ 排序库基因名 = %d" % len(hit))
print("  样例:", sorted(hit)[:25])
tf_tbl = pd.read_csv(os.path.join(DB, "motifs-v9-nr.flybase.tbl"), sep="\t", low_memory=False)
fly_tfs = set(tf_tbl["gene_name"].dropna().astype(str))
print("\n  我的数据中：能对上果蝇 TF 名 且 在排序库里 = %d 个 TF" % len(mysym & fly_tfs & feat))
pd.Series(sorted(hit)).to_csv(out + r"\scenic_gene_overlap.csv", index=False, header=False)
