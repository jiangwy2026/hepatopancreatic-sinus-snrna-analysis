import sys, io, os, pandas as pd
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"; DB = out + r"\scenic_db"
tf_tbl = pd.read_csv(os.path.join(DB, "motifs-v9-nr.flybase.tbl"), sep="\t", low_memory=False)
fly_tfs = set(tf_tbl["gene_name"].dropna().astype(str))
print("果蝇 TF 数:", len(fly_tfs), "| 样例:", sorted(fly_tfs)[:12])
feats = pd.read_csv(root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\1.Expression\expressions\blank\features.tsv.gz", sep="\t", header=None, names=["gid","sym","type"])
mysym = set(feats["sym"].astype(str))
hit = sorted(mysym & fly_tfs)
print("\n我的数据里能直接对上果蝇 TF 名的: %d 个" % len(hit))
print("  样例:", hit[:30])
pd.Series(hit).to_csv(out + r"\scenic_tf_overlap.csv", index=False, header=False)
print("\n(对照)我的数据中带 TF_family 注释的基因数: 见前次分析 263 个")
# 排名库的基因名空间
import pyarrow.feather as feather
f = os.path.join(DB, "dm6_v10_clust.regions_vs_motifs.rankings.feather")
try:
    with open(f, "rb") as fh:
        tbl = feather.read_table(fh, memory_map=True)
    print("\n=== 排序库 feather ===")
    print("  形状:", tbl.num_rows, "行 x", tbl.num_columns, "列")
    print("  行名前 5:", [str(x) for x in tbl.column(0)[:5]])
    print("  列名前 5:", [str(x) for x in tbl.column_names[:5]])
except Exception as e:
    print("\nfeather 读取(可能还在下载):", str(e)[:120])
