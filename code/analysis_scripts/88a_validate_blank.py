# -*- coding: utf-8 -*-
"""88a_validate_blank.py - check the custom velocity counts against the vendor filtered matrix."""
import os, gzip
import numpy as np
import scipy.sparse as sp
import scipy.io as sio

OUT = r"D:\桌面\single cell\analysis_output"
EXPR = os.path.join(r"D:\桌面\single cell", "ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result",
                    "ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result", "1.Expression", "expressions", "blank")
NG, NC = 26778, 11753
z = np.load(os.path.join(OUT, "velocity_counts_blank.npz"))
cells, genes, kinds = z["cells"].astype(np.int64), z["genes"].astype(np.int64), z["kinds"]
print("keys:", len(cells), "| spliced keys:", int((kinds == 0).sum()),
      "| unspliced:", int((kinds == 1).sum()), "| ambiguous:", int((kinds == 2).sum()))
def mat(k):
    sel = kinds == k
    return sp.coo_matrix((np.ones(int(sel.sum()), np.float32), (cells[sel], genes[sel])),
                         shape=(NC, NG)).tocsr()
S, U, A = mat(0), mat(1), mat(2)
T = (S + U + A).tocsr()
with gzip.open(os.path.join(EXPR, "features.tsv.gz"), "rt") as f:
    feat = [l.split("\t")[0] for l in f]
with gzip.open(os.path.join(EXPR, "barcodes.tsv.gz"), "rt") as f:
    bcs = [(l.strip()[:-2] if l.strip().endswith("-1") else l.strip()) for l in f if l.strip()]
vend = sio.mmread(os.path.join(EXPR, "matrix.mtx.gz")).tocsc()   # genes x cells
print("vendor matrix:", vend.shape, "| total UMI:", int(vend.sum()))
ours_tot = np.asarray(T.sum(axis=1)).ravel()
vend_tot = np.asarray(vend.sum(axis=0)).ravel()
ok = vend_tot > 0
print("cells with vendor UMI>0: %d/%d" % (ok.sum(), NC))
print("median UMI: vendor %.0f | ours(total) %.0f | ours(spliced) %.0f | ours(unspliced) %.0f" % (
    np.median(vend_tot[ok]), np.median(ours_tot[ok]),
    np.median(np.asarray(S.sum(axis=1)).ravel()[ok]), np.median(np.asarray(U.sum(axis=1)).ravel()[ok])))
print("pearson ours_total vs vendor: %.3f | spearman %.3f" % (
    np.corrcoef(ours_tot[ok], vend_tot[ok])[0, 1],
    sp.stats.spearmanr(ours_tot[ok], vend_tot[ok]).statistic if hasattr(sp, "stats") else float("nan")))
from scipy.stats import spearmanr
print("spearman total: %.3f | spearman spliced-only: %.3f" % (
    spearmanr(ours_tot[ok], vend_tot[ok]).statistic,
    spearmanr(np.asarray(S.sum(axis=1)).ravel()[ok], vend_tot[ok]).statistic))
print("median ratio ours/vendor: %.3f" % np.median(ours_tot[ok] / np.maximum(vend_tot[ok], 1)))
# per-gene check
vg = np.asarray(vend.sum(axis=1)).ravel()
og = np.asarray(T.sum(axis=0)).ravel()
osp = np.asarray(S.sum(axis=0)).ravel()
ogn = np.asarray(U.sum(axis=0)).ravel()
print("genes with vendor>0: %d | ours>0: %d | overlap %d" % ((vg > 0).sum(), (og > 0).sum(), ((vg > 0) & (og > 0)).sum()))
print("gene-level spearman: total %.3f | spliced %.3f" % (spearmanr(og, vg).statistic, spearmanr(osp, vg).statistic))
print("unspliced genes >0: %d | genes with unspliced share >10%%: %d" % (
    (ogn > 0).sum(), ((ogn / np.maximum(og, 1)) > 0.10).sum()))
order = np.argsort(-vg)[:12]
print("top vendor genes -> ours(total/spliced/unspliced):")
for i in order:
    print("   %-16s vendor %7d | ours %7d | spl %7d | uns %6d" % (feat[i], vg[i], og[i], osp[i], ogn[i]))
