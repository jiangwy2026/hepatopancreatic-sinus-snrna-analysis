# -*- coding: utf-8 -*-
"""
88_velocity_assemble.py - build spliced/unspliced/ambiguous layers from molecule keys and
attach them to the SoupX-corrected object; validate against the vendor matrices.
"""
import os, gzip, json
import numpy as np
import scipy.sparse as sp
import scipy.io as sio
import anndata as ad

ROOT = r"D:\桌面\single cell"
OUT = os.path.join(ROOT, "analysis_output")
EXPR = os.path.join(ROOT, "ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result",
                    "ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result", "1.Expression", "expressions")
NG = 26778
TOT = 19795
SAMPLES = [("blank", 0, 11753), ("treated", 11753, 8042)]
KINDS = ["spliced", "unspliced", "ambiguous"]

def read_barcodes(s):
    with gzip.open(os.path.join(EXPR, s, "barcodes.tsv.gz"), "rt") as f:
        return [(x[:-2] if x.endswith("-1") else x) for x in (l.strip() for l in f) if x]

# ---- molecule keys -> sparse matrices
mats = {k: np.zeros((TOT, NG), dtype=np.float32) for k in KINDS}
layers = {}
report = {}
for s, off, n in SAMPLES:
    z = np.load(os.path.join(OUT, "velocity_counts_%s.npz" % s))
    cells = z["cells"].astype(np.int64); genes = z["genes"].astype(np.int64); kinds = z["kinds"]
    per = {}
    for ki, kname in enumerate(KINDS):
        sel = kinds == ki
        coo = sp.coo_matrix((np.ones(int(sel.sum()), dtype=np.float32), (cells[sel], genes[sel])),
                            shape=(TOT, NG))
        c = coo.tocsr()
        layers.setdefault(kname, []).append(c)
        per[kname] = float(c.sum())
    report[s] = per
    print(s, "molecule counts:", {k: int(v) for k, v in per.items()},
          "| unspliced share %.1f%%" % (100 * per["unspliced"] / max(sum(per.values()), 1)))
S = sum(layers["spliced"]).tocsr()[:TOT]
U = sum(layers["unspliced"]).tocsr()[:TOT]
A = sum(layers["ambiguous"]).tocsr()[:TOT]
TOTM = (S + U + A).tocsr()
print("total molecules:", int(TOTM.sum()), "| genes with unspliced>0:",
      int((np.asarray(U.sum(axis=0)).ravel() > 0).sum()),
      "| genes with spliced>0:", int((np.asarray(S.sum(axis=0)).ravel() > 0).sum()))

# ---- validation against vendor matrices
for s, off, n in SAMPLES:
    m = sio.mmread(os.path.join(EXPR, s, "matrix.mtx.gz"))   # genes x cells
    vend = np.asarray(m.sum(axis=0)).ravel().astype(np.float64)
    ours = np.asarray(TOTM[off:off + n].sum(axis=1)).ravel()
    ours_sp = np.asarray(S[off:off + n].sum(axis=1)).ravel()
    ours_un = np.asarray(U[off:off + n].sum(axis=1)).ravel()
    ok = vend > 0
    r_all = np.corrcoef(ours[ok], vend[ok])[0, 1]
    r_sp = np.corrcoef(ours_sp[ok], vend[ok])[0, 1]
    ratio = np.median(ours[ok] / np.maximum(vend[ok], 1))
    print("%s: vendor median UMI %.0f vs ours %.0f | pearson(all) %.3f | pearson(exonic-only) %.3f | median ratio %.3f" % (
        s, np.median(vend[ok]), np.median(ours[ok]), r_all, r_sp, ratio))
    print("   unspliced per cell: median %.0f | mean %.0f | cells with >0: %.1f%%" % (
        np.median(ours_un[ok]), ours_un[ok].mean(), 100 * (ours_un[ok] > 0).mean()))
    report[s].update({"vendor_median_umi": float(np.median(vend[ok])), "ours_median_total": float(np.median(ours[ok])),
                      "pearson_all": float(r_all), "pearson_exonic": float(r_sp), "median_ratio": float(ratio),
                      "median_unspliced_per_cell": float(np.median(ours_un[ok])),
                      "pct_cells_with_unspliced": float(100 * (ours_un[ok] > 0).mean())})

# ---- attach to object
a = ad.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))
lookup = {}
for s, off, n in SAMPLES:
    for i, b in enumerate(read_barcodes(s)):
        lookup[s + "_" + b] = off + i
idx = np.array([lookup.get(x, -1) for x in a.obs_names], dtype=np.int64)
print("cells mapped to BAM:", int((idx >= 0).sum()), "of", a.n_obs)
assert (idx >= 0).all()
a.layers["spliced"] = S[idx]
a.layers["unspliced"] = U[idx]
a.layers["ambiguous"] = A[idx]
a.uns["velocity_source"] = "custom BGZF/struct counter over CellRanger 8 possorted BAM, UMI-deduplicated, GTF-based exon/intron classification"
print("obs groups:")
if "group" in a.obs.columns:
    print(a.obs["group"].value_counts().to_string())
sp2 = np.asarray(a.layers["spliced"].sum(axis=1)).ravel()
un2 = np.asarray(a.layers["unspliced"].sum(axis=1)).ravel()
tot2 = sp2 + un2 + np.asarray(a.layers["ambiguous"].sum(axis=1)).ravel()
print("object total molecules median %.0f | spliced median %.0f | unspliced median %.0f" % (
    np.median(tot2), np.median(sp2), np.median(un2)))
a.write_h5ad(os.path.join(OUT, "adata_velocity_raw.h5ad"))
with open(os.path.join(OUT, "VELOCITY_COUNT_REPORT.json"), "w") as f:
    json.dump(report, f, indent=1)
print("saved adata_velocity_raw.h5ad")
