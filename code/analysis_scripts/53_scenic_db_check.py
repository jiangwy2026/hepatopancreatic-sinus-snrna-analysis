import sys, io, os, pandas as pd, numpy as np
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"; DB = out + r"\scenic_db"
import pyarrow.feather as feather
f = os.path.join(DB, "dm6_v10_clust.regions_vs_motifs.rankings.feather")
tbl = feather.read_table(f, memory_map=True)
print("=== 排序库 ===")
print("  形状: %d 行 x %d 列" % (tbl.num_rows, tbl.num_columns))
print("  行名(基因) 前 8:", [str(x) for x in tbl.column(0)[:8]])
print("  列名(motif) 前 5:", [str(x) for x in tbl.column_names[:5]])
feat = set(str(x) for x in tbl.column(0))
print("  基因总数: %d" % len(feat))
# 我的数据里有多少基因能对上
feats = pd.read_csv(root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\1.Expression\expressions\blank\features.tsv.gz", sep="\t", header=None, names=["gid","sym","type"])
mysym = set(feats["sym"].astype(str))
hit = mysym & feat
print("\n  我的基因符号与排序库基因名的交集: %d" % len(hit))
print("  样例:", sorted(hit)[:25])
pd.Series(sorted(hit)).to_csv(out + r"\scenic_gene_overlap.csv", index=False, header=False)
