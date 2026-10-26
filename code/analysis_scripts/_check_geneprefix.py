# -*- coding: utf-8 -*-
import os, gzip, re, numpy as np, scipy.io as sio
EXPR = os.path.join(r"D:\桌面\single cell", "ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result",
                    "ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result", "1.Expression", "expressions", "blank")
GTF = r"E:\workspace\genome2\genome\GENOME_GTF\GCF_003789085.1.gtf.gz"
with gzip.open(os.path.join(EXPR, "features.tsv.gz"), "rt") as f:
    rows = [l.rstrip("\n").split("\t") for l in f]
gene_pref = [r for r in rows if r[0].startswith("gene-")]
print("features starting with 'gene-':", len(gene_pref), "examples:", [r[0] for r in gene_pref[:3]])
gtf_ids = set()
with gzip.open(GTF, "rt") as f:
    for line in f:
        if line[0] != "#" and "\tgene\t" in line:
            m = re.search(r'gene_id "([^"]+)"', line)
            if m: gtf_ids.add(m.group(1))
print("GTF gene ids:", len(gtf_ids))
stripped = [r for r in gene_pref if r[0][5:] in gtf_ids]
print("'gene-' features whose stripped id IS in GTF:", len(stripped))
vend = sio.mmread(os.path.join(EXPR, "matrix.mtx.gz"))
vg = np.asarray(vend.sum(axis=1)).ravel()
idx = {r[0]: i for i, r in enumerate(rows)}
mass_all = sum(vg[idx[r[0]]] for r in gene_pref)
mass_fix = sum(vg[idx[r[0]]] for r in stripped)
print("vendor UMI mass in 'gene-' features: %d | recoverable by stripping: %d" % (mass_all, mass_fix))
# how many matrix features fail all our mapping rules?
n_map = 0
for r in rows:
    c0, c1 = r[0], r[1]
    ok = False
    if c0.startswith("ncbi_"):
        cand = "LOC" + c0[5:]
        ok = cand in gtf_ids
    if not ok: ok = c1 in gtf_ids
    if not ok: ok = c0 in gtf_ids
    if not ok and c0.startswith("gene-"):
        ok = c0[5:] in gtf_ids
    if ok: n_map += 1
print("total features:", len(rows), "| mappable with current+strip rules:", n_map)
