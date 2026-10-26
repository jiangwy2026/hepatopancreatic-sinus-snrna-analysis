import sys, io, os, glob
import numpy as np, pandas as pd, scanpy as sc
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 250); sc.settings.verbosity = 0
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"; GEO = out + r"\geo_ref"
files = sorted(glob.glob(os.path.join(GEO, "GSM*_filtered_feature_bc_matrix.h5")))
print("=== GSE306674 各样本 ===")
tot = {}
A0 = None
for f in files:
    name = os.path.basename(f).split("_")[1]
    a = sc.read_10x_h5(f)
    print("  %-4s %6d 细胞 x %6d 基因  | 中位 counts %.0f | 中位基因数 %.0f | 首基因 %s" %
          (name, a.n_obs, a.n_vars, np.median(np.asarray(a.X.sum(1)).ravel()), np.median(np.asarray((a.X>0).sum(1)).ravel()), a.var_names[0]), flush=True)
    tot[name] = a
    if A0 is None: A0 = a
print("\n总细胞数:", sum(a.n_obs for a in tot.values()))
print("\n=== 基因注释格式对比 ===")
print("  GSE306674 前 5 个基因:", list(A0.var_names[:5]))
mine = pd.read_csv(root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\1.Expression\expressions\blank\features.tsv.gz", sep="\t", header=None)[0].astype(str).tolist()
print("  你的数据  前 5 个基因:", mine[:5])
ref = set(A0.var_names); my = set(mine)
print("\n  参考基因 %d | 你的基因 %d | 交集 %d (占你数据 %.1f%%)" % (len(ref), len(my), len(ref & my), len(ref & my)/len(my)*100))
print("  GSE306674 var 列:", list(A0.var.columns), "| 你的 var 列: ['gid/sym (来自 features.tsv)']")
