import sys, io, os, gzip
import numpy as np, pandas as pd, scipy.sparse as sp
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 260)
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
FIL = root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\1.Expression\expressions"
def rb(p):
    with gzip.open(p, "rt") as f: return [l.rstrip("\n") for l in f]
feats = pd.read_csv(FIL + r"\blank\features.tsv.gz", sep="\t", header=None, names=["gid","sym","type"])
sym = feats["sym"].astype(str).values; G = len(sym)
sym2idx = {}
for i, x in enumerate(sym): sym2idx.setdefault(x, i)
mk = pd.read_csv(out + r"\markers_vendor_scanpy.csv"); mk["gene"] = mk["symbol"].fillna(mk["names"]).astype(str)
top = mk[(mk.pvals_adj < 0.01) & (mk.logfoldchanges > 0.5) & (mk.pct_nz_group > 0.20)]
mask_gene = {str(c): np.array([sym2idx[g] for g in set(top[top.cluster.astype(str)==str(c)].gene) if g in sym2idx], dtype=int)
             for c in sorted(top.cluster.unique(), key=lambda x:int(x))}
print("各群用于屏蔽的基因数:", {k: len(v) for k,v in mask_gene.items()})

cl = pd.read_csv(root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\2.Cluster\2.cluster\Cells.cluster.list.xls", sep="\t")
cl["Cells"]=cl["Cells"].astype(str); cl["Cluster"]=cl["Cluster"].astype(str)
summary = []
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
    a = np.zeros(F.shape[0])
    for c in sorted(set(clu), key=int):
        idx = clu[clu==c].index.to_numpy()
        keep = np.ones(G, dtype=bool); keep[mask_gene.get(str(c), np.array([], dtype=int))] = False
        sv = soup[keep]; sv = sv - sv.mean()
        den = (sv*sv).sum()
        # 对该群细胞：a_i = (X_keep · sv) / ||sv||^2
        blk = F[idx][:, keep]
        a[idx] = np.asarray(blk @ sv).ravel() / den
    rho = np.clip(a/np.maximum(tot,1), 0, 0.95)
    print("  rho 中位=%.3f | 25%%=%.3f | 75%%=%.3f | a 中位=%.0f (细胞总 UMI 中位 %.0f)" %
          (np.median(rho), np.percentile(rho,25), np.percentile(rho,75), np.median(a), np.median(tot)))
    print("  按群 rho 中位:", {c: round(float(np.median(rho[clu[clu==c].index.to_numpy()])),3) for c in sorted(set(clu), key=int)})
    Cf = np.zeros(F.shape, dtype=np.float32); soup32 = soup.astype(np.float32); a32 = (rho*tot).astype(np.float32)
    for st in range(0, F.shape[0], 2000):
        en = min(st+2000, F.shape[0])
        Cf[st:en] = np.maximum(F[st:en].toarray().astype(np.float32) - a32[st:en,None]*soup32[None,:], 0)
    Cf = sp.csr_matrix(Cf); sp.save_npz(out + r"\_tmp_corrected_%s.npz" % s, Cf); np.save(out + r"\_tmp_rho_%s.npy" % s, rho)
    print("  总 counts %.0f -> %.0f (减少 %.1f%%)" % (F.sum(), Cf.sum(), (1-Cf.sum()/F.sum())*100))
    summary.append(dict(sample=s, rho_median=float(np.median(rho)), rho_q25=float(np.percentile(rho,25)),
                        rho_q75=float(np.percentile(rho,75)), counts_before=float(F.sum()), counts_after=float(Cf.sum())))
pd.DataFrame(summary).to_csv(out + r"\soup_correction_summary.csv", index=False)
print("\nwrote soup_correction_summary.csv")
