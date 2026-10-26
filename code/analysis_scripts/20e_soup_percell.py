import sys, io, os, gzip
import numpy as np, pandas as pd, scipy.sparse as sp
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 260)
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
FIL = root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\1.Expression\expressions"
def rb(p):
    with gzip.open(p, "rt") as f: return [l.rstrip("\n") for l in f]
feats = pd.read_csv(FIL + r"\blank\features.tsv.gz", sep="\t", header=None, names=["gid","sym","type"])
sym = feats["sym"].astype(str).values
sym2idx = {}
for i, x in enumerate(sym): sym2idx.setdefault(x, i)

# 高特异性 marker：只在一个群里表达，在其余群里几乎不表达
mk = pd.read_csv(out + r"\markers_vendor_scanpy.csv"); mk["gene"] = mk["symbol"].fillna(mk["names"]).astype(str)
strict = mk[(mk.pvals_adj < 0.01) & (mk.logfoldchanges > 1.5) & (mk.pct_nz_group > 0.40) & (mk.pct_nz_reference < 0.30)]
print("高特异性 marker 数（按群）:", strict.groupby(strict.cluster.astype(str)).size().to_dict())
mset = {str(c): sorted(set(strict[strict.cluster.astype(str)==str(c)].gene)) for c in sorted(strict.cluster.unique(), key=lambda x:int(x))}
for k, v in mset.items():
    print("   C%-2s %d 个" % (k, len(v)))

cl = pd.read_csv(root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\2.Cluster\2.cluster\Cells.cluster.list.xls", sep="\t")
cl["Cells"]=cl["Cells"].astype(str); cl["Cluster"]=cl["Cluster"].astype(str)
report = []
for s in ("blank","treated"):
    print("\n################ %s ################" % s)
    F = sp.load_npz(out + r"\_tmp_filt_%s.npz" % s).tocsr()
    soup = np.load(out + r"\_tmp_soup_%s.npy" % s)
    tot = np.asarray(F.sum(1)).ravel()
    fb = [b.rsplit("-",1)[0] for b in rb(FIL + "\\" + s + r"\barcodes.tsv.gz")]
    key = pd.Series(range(len(fb)), index=[s + "_" + b for b in fb])
    cc = cl[cl.Cells.str.startswith(s + "_")].copy(); cc["idx"] = key.reindex(cc.Cells).values
    cc = cc.dropna(subset=["idx"]); cc["idx"] = cc["idx"].astype(int)
    clu = cc.set_index("idx")["Cluster"]
    # 逐细胞：用“其他群的高特异 marker”估计环境总量 a_i
    a = np.zeros(F.shape[0])
    for c in sorted(set(clu), key=int):
        idx = clu[clu==c].index.to_numpy()
        gs = [sym2idx[x] for k, v in mset.items() if k != str(c) for x in v if x in sym2idx]
        if not gs: continue
        gs = np.array(gs)
        num = np.asarray(F[idx][:, gs].sum(1)).ravel()
        den = soup[gs].sum()
        a[idx] = num / max(den, 1e-12)
    rho = a / np.maximum(tot, 1)
    rho = np.clip(rho, 0, 0.95)
    print("  环境总量 a_i 中位=%.0f | 污染率 rho 中位=%.3f (25%%=%.3f, 75%%=%.3f)" % (np.median(a), np.median(rho), np.percentile(rho,25), np.percentile(rho,75)))
    print("  按群 rho 中位:", {c: round(float(np.median(rho[clu[clu==c].index.to_numpy()])),3) for c in sorted(set(clu), key=int)})
    # 校正
    a_clip = rho * tot
    C = F.astype(np.float64).tolil()
    step = 2000
    Cf = np.zeros(F.shape, dtype=np.float32)
    for st in range(0, F.shape[0], step):
        en = min(st+step, F.shape[0])
        blk = F[st:en].toarray().astype(np.float32)
        blk = np.maximum(blk - (a_clip[st:en][:,None]*soup[None,:]).astype(np.float32), 0)
        Cf[st:en] = blk
    Cf = sp.csr_matrix(Cf)
    sp.save_npz(out + r"\_tmp_corrected_%s.npz" % s, Cf)
    np.save(out + r"\_tmp_rho_%s.npy" % s, rho)
    print("  总 counts %.0f -> %.0f (减少 %.1f%%) | 非零元素 %.1fM -> %.1fM" % (F.sum(), Cf.sum(), (1-Cf.sum()/F.sum())*100, F.nnz/1e6, Cf.nnz/1e6))
    report.append(dict(sample=s, rho_median=float(np.median(rho)), a_median=float(np.median(a)), counts_before=float(F.sum()), counts_after=float(Cf.sum())))
pd.DataFrame(report).to_csv(out + r"\soup_correction_summary.csv", index=False)
print("\nwrote soup_correction_summary.csv")
