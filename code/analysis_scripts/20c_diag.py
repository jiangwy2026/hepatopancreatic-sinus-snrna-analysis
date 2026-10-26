import sys, io, os, gzip
import numpy as np, pandas as pd, scipy.sparse as sp
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 250)
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
FIL = root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\1.Expression\expressions"
def rb(p):
    with gzip.open(p, "rt") as f: return [l.rstrip("\n") for l in f]
feats = pd.read_csv(FIL + r"\blank\features.tsv.gz", sep="\t", header=None, names=["gid","sym","type"])
sym = feats["sym"].astype(str).values; G = len(sym)
soup = np.load(out + r"\_tmp_soup_blank.npy")
o = np.argsort(-soup)
print("=== soup 谱集中度 ===")
for k in (1, 10, 50, 100, 500):
    print("  top %-4d 基因占 soup 的 %.1f%%" % (k, soup[o[:k]].sum()*100))
print("  top10:", ", ".join(sym[o[:10]]))
print("  soup 非零基因数: %d (%.0f%%)" % ((soup>0).sum(), (soup>0).mean()*100))

F = sp.load_npz(out + r"\_tmp_filt_blank.npz" % ()).tocsr()
cl = pd.read_csv(root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\2.Cluster\2.cluster\Cells.cluster.list.xls", sep="\t")
cl["Cells"]=cl["Cells"].astype(str); cl["Cluster"]=cl["Cluster"].astype(str)
fb = [b.rsplit("-",1)[0] for b in rb(FIL + r"\blank\barcodes.tsv.gz")]
key = pd.Series(range(len(fb)), index=["blank_" + b for b in fb])
cc = cl[cl.Cells.str.startswith("blank_")].copy(); cc["idx"] = key.reindex(cc.Cells).values
cc = cc.dropna(subset=["idx"]); cc["idx"] = cc["idx"].astype(int)
clu = cc.set_index("idx")["Cluster"]
tot = np.asarray(F.sum(1)).ravel()

mk = pd.read_csv(out + r"\markers_vendor_scanpy.csv"); mk["gene"] = mk["symbol"].fillna(mk["names"]).astype(str)
sig = mk[(mk.pvals_adj < 0.01) & (mk.logfoldchanges > 1.0) & (mk.pct_nz_group > 0.30)]
markers = {str(c): set(sig[sig.cluster.astype(str)==str(c)].gene) for c in sorted(sig.cluster.unique(), key=lambda x:int(x))}
sym2idx = {}
for i, x in enumerate(sym): sym2idx.setdefault(x, i)

print("\n=== rho 估计诊断（blank, 亚群 0 / 1 / 4）===")
for c in ["0","1","4"]:
    idx = clu[clu==c].index.to_numpy(); U = tot[idx].mean()
    sub = F[idx]
    det = np.asarray((sub>0).mean(0)).ravel()             # 检出率
    obs = np.asarray(sub.mean(0)).ravel().astype(float)   # 平均 counts
    other = set().union(*[v for k,v in markers.items() if k!=c]) - markers[c]
    g_all = np.array([sym2idx[x] for x in other if x in sym2idx])
    g_strict = g_all[det[g_all] < 0.05]                   # 严格“不表达”
    print("\n  --- 亚群 %s | 细胞 %d | 平均 UMI %.0f" % (c, len(idx), U))
    print("      其他群 marker 基因 %d 个；其中在本群检出率<5%% 的 %d 个" % (len(g_all), len(g_strict)))
    for tag, g in (("全部", g_all), ("严格", g_strict)):
        if len(g) < 10: print("      [%s] 基因太少，跳过" % tag); continue
        a = soup[g]*U; ratio = obs[g]/np.maximum(a, 1e-12)
        ls = float((obs[g]*a).sum()/(a*a).sum())
        print("      [%s] 最小二乘 rho=%.3f | 比值中位数=%.3f | 25分位=%.3f | 75分位=%.3f" %
              (tag, ls, np.median(ratio), np.percentile(ratio,25), np.percentile(ratio,75)))
    top = np.argsort(-obs[g_strict])[:8] if len(g_strict)>=8 else np.argsort(-obs[g_all])[:8]
    for t in top[:6]:
        i = (g_strict if len(g_strict)>=8 else g_all)[t]
        print("        %-18s obs=%7.2f  预期soup=%7.2f  比值=%.2f" % (sym[i], obs[i], soup[i]*U, obs[i]/max(soup[i]*U,1e-12)))
