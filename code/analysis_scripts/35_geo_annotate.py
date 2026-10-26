import sys, io, os, re
import numpy as np, pandas as pd
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 260)
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
feats = pd.read_csv(root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\1.Expression\expressions\blank\features.tsv.gz", sep="\t", header=None, names=["gid","sym","type"])
num = lambda s: re.sub(r"^(ncbi_|LOC|gene-)", "", str(s))
id2sym = {num(g): s for g, s in zip(feats["gid"], feats["sym"])}
mk = pd.read_csv(out + r"\geo_ref_cluster_markers.csv")
mk["sym"] = mk["names"].astype(str).map(lambda x: id2sym.get(num(x), "-"))
mk.to_csv(out + r"\geo_ref_cluster_markers_sym.csv", index=False)
print("=== 参考数据各群 top12 marker（带符号）===")
for g in sorted(mk["cl"].unique(), key=lambda x: int(x)):
    d = mk[mk["cl"] == g].sort_values("scores", ascending=False).head(12)
    hit = sum(1 for s in d["sym"] if s != "-")
    print("  R%-3s: %s   [匹配到符号 %d/12]" % (g, ", ".join(d["sym"].astype(str).tolist()), hit))
# 关键 marker 在各参考群的表达
probe = {"Prss1":"F细胞(胰蛋白酶)","Try3":"F细胞","Cpa2":"F细胞","LCP2":"F细胞",
         "ACE":"R细胞","Oacyl":"R细胞","LIPF":"R细胞","Mal-A2":"B细胞","MALRD1":"B细胞",
         "A2M":"颗粒血细胞","PPAF2":"血细胞","LOC113823945":"增殖血细胞","ccnd2":"增殖",
         "LOC113809538":"表皮","RASGRF1":"神经","PEN-3":"抗菌肽"}
sym2num = {v: k for k, v in id2sym.items()}
print("\n=== 关键 marker 在参考各群的表达位次（由 marker 表推断）===")
for s, tag in probe.items():
    n = sym2num.get(s)
    if n is None: print("  %-14s (%s) 未在参考基因集中" % (s, tag)); continue
    rows = mk[mk["names"].astype(str).map(num) == n]
    if rows.empty: print("  %-14s (%s) 未进入任何群的 top marker" % (s, tag)); continue
    rows = rows.sort_values("scores", ascending=False).head(3)
    print("  %-14s (%s) -> %s" % (s, tag, "; ".join("R%s score=%.1f" % (r["cl"], r["scores"]) for _, r in rows.iterrows())))
