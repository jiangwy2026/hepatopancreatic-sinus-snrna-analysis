import sys, io, os, gzip, time
import numpy as np, pandas as pd, scipy.sparse as sp
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
RAW = root + r"\ADRL26080060_raw"
FIL = root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\1.Expression\expressions"

def read_barcodes(p):
    with gzip.open(p, "rt") as f: return [l.rstrip("\n") for l in f]

def load_mtx(p, n_gene, n_cell, dtype=np.float32):
    t0 = time.time(); rows, cols = [], []
    dfs = []
    for ch in pd.read_csv(p, sep=" ", skiprows=3, header=None, names=["g","c","x"],
                          dtype={"g":np.int32,"c":np.int32,"x":np.int32},
                          compression="gzip", engine="c", chunksize=6_000_000):
        dfs.append(ch)
    d = pd.concat(dfs, ignore_index=True); del dfs
    M = sp.csr_matrix((d["x"].to_numpy(dtype=dtype), (d["c"].to_numpy()-1, d["g"].to_numpy()-1)),
                      shape=(n_cell, n_gene))
    print("    loaded %s: %d x %d nnz=%d  (%.0fs)" % (os.path.basename(os.path.dirname(p)), M.shape[0], M.shape[1], M.nnz, time.time()-t0))
    return M

feats = pd.read_csv(FIL + r"\blank\features.tsv.gz", sep="\t", header=None, names=["gid","sym","type"])
NG = len(feats)
res = {}
for s in ("blank","treated"):
    print("=== 样本 %s" % s)
    rb = read_barcodes(RAW + "\\" + s + r"\raw_feature_bc_matrix\barcodes.tsv.gz")
    fb = read_barcodes(FIL + "\\" + s + r"\barcodes.tsv.gz")
    R = load_mtx(RAW + "\\" + s + r"\raw_feature_bc_matrix\matrix.mtx.gz", NG, len(rb))
    F = load_mtx(FIL + "\\" + s + r"\matrix.mtx.gz", NG, len(fb))
    tot_all = np.asarray(R.sum(1)).ravel()
    tot_f = np.asarray(F.sum(1)).ravel()
    print("    raw barcodes=%d | 过滤后 barcodes=%d" % (R.shape[0], F.shape[0]))
    print("    raw 总 UMI 分位: 50%%=%.0f 90%%=%.0f 99%%=%.0f max=%.0f" % (np.percentile(tot_all,50), np.percentile(tot_all,90), np.percentile(tot_all,99), tot_all.max()))
    print("    过滤矩阵总 UMI 分位: 中位=%.0f max=%.0f" % (np.median(tot_f), tot_f.max()))
    sp.arange = None
    res[s] = dict(total_all=tot_all, total_f=tot_f, n_gene=NG, raw_path=RAW+"\\"+s, n_raw=len(rb))
    sp.save_npz(out + r"\_tmp_raw_%s.npz" % s, R.tocsr())
    sp.save_npz(out + r"\_tmp_filt_%s.npz" % s, F.tocsr())
    np.save(out + r"\_tmp_totall_%s.npy" % s, tot_all)
    print("    saved temp matrices")
print("\nDONE")
