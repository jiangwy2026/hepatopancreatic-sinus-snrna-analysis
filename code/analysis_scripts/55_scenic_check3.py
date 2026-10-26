import sys, io, os, pandas as pd, numpy as np
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"; DB = out + r"\scenic_db"
import pyarrow.feather as feather
tbl = feather.read_table(os.path.join(DB, "dm6_v10_clust.genes_vs_motifs.rankings.feather"), memory_map=True)
names = [str(x) for x in tbl.column_names]
print("=== 结构 ===")
print("  行数(motif) = %d | 列数(基因) = %d" % (tbl.num_rows, tbl.num_columns))
print("  列名(基因) 前 8:", names[:8])
feat = set(names)
feats = pd.read_csv(root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\1.Expression\expressions\blank\features.tsv.gz", sep="\t", header=None, names=["gid","sym","type"])
mysym = set(feats["sym"].astype(str))
hit = mysym & feat
print("\n  排序库基因数: %d" % len(feat))
print("  我的基因符号 ∩ 排序库 = %d 个 (占我全部基因的 %.1f%%)" % (len(hit), len(hit)/len(mysym)*100))
print("  样例:", sorted(hit)[:20])
tf_tbl = pd.read_csv(os.path.join(DB, "motifs-v9-nr.flybase.tbl"), sep="\t", low_memory=False)
fly_tfs = set(tf_tbl["gene_name"].dropna().astype(str))
print("\n  果蝇 TF 数: %d | 其中出现在我数据里的: %d" % (len(fly_tfs), len(mysym & fly_tfs)))
print("  且同时在排序库里的: %d" % len(mysym & fly_tfs & feat))
pd.Series(sorted(hit)).to_csv(out + r"\scenic_gene_overlap.csv", index=False, header=False)
pd.Series(sorted(mysym & fly_tfs & feat)).to_csv(out + r"\scenic_TF_list.csv", index=False, header=False)
