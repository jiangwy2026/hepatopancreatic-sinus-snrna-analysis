import sys, io, os, gzip
import numpy as np, pandas as pd, scipy.sparse as sp
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 250)
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
FIL = root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\1.Expression\expressions"

def rb(p):
    with gzip.open(p, "rt") as f: return [l.rstrip("\n") for l in f]

feats = pd.read_csv(FIL + r"\blank\features.tsv.gz", sep="\t", header=None, names=["gid","sym","type"])
sym = feats["sym"].astype(str).values
cl = pd.read_csv(root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\2.Cluster\2.cluster\Cells.cluster.list.xls", sep="\t")
cl["Cells"] = cl["Cells"].astype(str); cl["Cluster"] = cl["Cluster"].astype(str)
mk = pd.read_csv(out + r"\markers_vendor_scanpy.csv")
mk["gene"] = mk["symbol"].fillna(mk["names"]).astype(str)
sig = mk[(mk.pvals_adj < 0.01) & (mk.logfoldchanges > 1.0) & (mk.pct_nz_group > 0.30)]
markers = {str(c): set(sig[sig.cluster.astype(str) == str(c)].gene) for c in sorted(sig.cluster.unique(), key=lambda x: int(x))}
print("各群标记基因数:", {k: len(v) for k, v in markers.items()})
sym2idx = {}
for i, x in enumerate(sym): sym2idx.setdefault(x, i)

summary = []
for s in ("blank","treated"):
    print("\n########## %s ##########" % s)
    F = sp.load_npz(out + r"\_tmp_filt_%s.npz" % s).tocsr()
    soup = np.load(out + r"\_tmp_soup_%s.npy" % s)
    tot_f = np.asarray(F.sum(1)).ravel()
    fb = [b.rsplit("-",1)[0] for b in rb(FIL + "\\" + s + r"\barcodes.tsv.gz")]
    key = pd.Series(range(len(fb)), index=[s + "_" + b for b in fb])
    cc = cl[cl.Cells.str.startswith(s + "_")].copy()
    cc["idx"] = key.reindex(cc.Cells).values
    cc = cc.dropna(subset=["idx"]); cc["idx"] = cc["idx"].astype(int)
    clu = cc.set_index("idx")["Cluster"]
    print("  匹配到分群的细胞: %d" % len(clu))
    rows = []
    for c in markers:
        idx = clu[clu == c].index.to_numpy()
        if len(idx) < 30: continue
        U_c = tot_f[idx].mean()
        obs = np.asarray(F[idx].mean(0)).ravel().astype(np.float64)
        other = (set().union(*[v for k, v in markers.items() if k != c]) - markers[c])
        g = np.array([sym2idx[x] for x in other if x in sym2idx], dtype=int)
        if len(g) < 20: continue
        a = soup[g] * U_c
        rho = float((obs[g] * a).sum() / max((a * a).sum(), 1e-12))
        rho = min(max(rho, 0.0), 0.9)
        rows.append(dict(sample=s, cluster=c, n_cells=len(idx), UMI_mean=round(U_c,1), rho=round(rho,4), n_marker_used=len(g)))
    d = pd.DataFrame(rows); summary.append(d)
    print(d.to_string(index=False))
    print("  平均 rho = %.3f | 中位 rho = %.3f" % (d.rho.mean(), d.rho.median()))
    C = F.astype(np.float64)
    for _, r in d.iterrows():
        idx = clu[clu == str(r["cluster"])].index.to_numpy()
        factor = (float(r["rho"]) * tot_f[idx])[:, None] * soup[None, :]
        C[idx] = np.maximum(C[idx] - factor, 0.0)
    sp.save_npz(out + r"\_tmp_corrected_%s.npz" % s, sp.csr_matrix(C.astype(np.float32)))
    print("  总 counts %.0f -> %.0f (减少 %.1f%%)" % (F.sum(), C.sum(), (1 - C.sum()/F.sum())*100))
pd.concat(summary).to_csv(out + r"\soupx_rho_per_cluster.csv", index=False)
print("\nwrote soupx_rho_per_cluster.csv")
