# -*- coding: utf-8 -*-
import os, gzip, re, numpy as np, scipy.sparse as sp, scipy.io as sio
from array import array
OUT = r"D:\桌面\single cell\analysis_output"
GTF = r"E:\workspace\genome2\genome\GENOME_GTF\GCF_003789085.1.gtf.gz"
EXPR = os.path.join(r"D:\桌面\single cell", "ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result",
                    "ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result", "1.Expression", "expressions", "blank")
# gene coordinates of the suspicious genes
loc = {}
with gzip.open(GTF, "rt") as f:
    for line in f:
        if line[0] == "#": continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] != "gene": continue
        m = re.search(r'gene_id "([^"]+)"', p[8])
        if m and m.group(1) in ("LOC113822661", "LOC113803873", "LOC113808997"):
            loc[m.group(1)] = (p[0], int(p[3]) - 1, int(p[4]))
print("gene loci:", loc)
# all genes on those contigs
byc = {}
with gzip.open(GTF, "rt") as f:
    for line in f:
        if line[0] == "#": continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] != "gene": continue
        m = re.search(r'gene_id "([^"]+)"', p[8])
        if not m: continue
        if p[0] in {v[0] for v in loc.values()}:
            byc.setdefault(p[0], []).append((int(p[3]) - 1, int(p[4]), m.group(1)))
for c, rows in byc.items():
    rows.sort()
    print("contig %s has %d genes; total span %d..%d" % (c, len(rows), rows[0][0], rows[-1][1]))
    for g in loc:
        if loc[g][0] == c:
            s, e = loc[g][1], loc[g][2]
            ov = [r for r in rows if r[1] > s and r[0] < e]
            print("   %s (%d-%d) overlaps %d genes: %s" % (g, s, e, len(ov), [(r[2], r[0], r[1]) for r in ov][:6]))
# genes with big disagreement
z = np.load(os.path.join(OUT, "velocity_counts_blank.npz"))
cells, genes, kinds = z["cells"].astype(np.int64), z["genes"].astype(np.int64), z["kinds"]
sel = kinds <= 1
T = sp.coo_matrix((np.ones(int(sel.sum()), np.float32), (cells[sel], genes[sel])), shape=(11753, 26778)).tocsr()
og = np.asarray(T.sum(axis=0)).ravel()
vend = sio.mmread(os.path.join(EXPR, "matrix.mtx.gz"))
vg = np.asarray(vend.sum(axis=1)).ravel()
with gzip.open(os.path.join(EXPR, "features.tsv.gz"), "rt") as f:
    feat = [l.rstrip("\n").split("\t") for l in f]
ratio = np.where(vg >= 100, og / np.maximum(vg, 1), np.nan)
bad = np.where(np.isfinite(ratio) & (ratio < 0.5))[0]
print("\ngenes with vendor>=100 UMI: %d | ratio<0.5: %d | ratio>2: %d" % (
    int((vg >= 100).sum()), len(bad), int((np.isfinite(ratio) & (ratio > 2)).sum())))
bad = bad[np.argsort(-vg[bad])][:15]
print("worst-covered highly expressed genes (vendor / ours / ratio):")
for i in bad:
    print("   %-16s %-14s vendor %7d ours %7d ratio %.2f" % (feat[i][0], feat[i][1], vg[i], og[i], ratio[i]))
