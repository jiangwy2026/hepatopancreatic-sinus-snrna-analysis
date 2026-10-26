import sys, io, os, pandas as pd, numpy as np
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"; DB = out + r"\scenic_db"
# 1) motif 注释表里的果蝇 TF 列表
tf_tbl = pd.read_csv(os.path.join(DB, "motifs-v9-nr.flybase.tbl"), sep="\t")
print("motif 注释表:", tf_tbl.shape, "| 列:", list(tf_tbl.columns))
col_tf = [c for c in tf_tbl.columns if c.lower() in ("tf","gene_name","#motif_id","motif_id")]
print("前 3 行:"); print(tf_tbl.head(3).to_string())
fly_tfs = set(tf_tbl["TF"].astype(str)) if "TF" in tf_tbl.columns else set()
print("\n果蝇 TF 总数:", len(fly_tfs))
print("样例:", sorted(list(fly_tfs))[:15])
# 2) 我的基因符号
feats = pd.read_csv(root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\1.Expression\expressions\blank\features.tsv.gz", sep="\t", header=None, names=["gid","sym","type"])
mysym = set(feats["sym"].astype(str))
print("\n我的基因符号数:", len(mysym))
hit = mysym & fly_tfs
print("与果蝇 TF 名直接匹配的基因数: %d" % len(hit))
print("  样例:", sorted(list(hit))[:25])
pd.Series(sorted(hit)).to_csv(out + r"\scenic_tf_overlap.csv", index=False, header=False)
