# -*- coding: utf-8 -*-
import os, gzip, re, collections
OUT = r"D:\桌面\single cell\analysis_output"
EXPR = os.path.join(r"D:\桌面\single cell", "ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result",
                    "ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result", "1.Expression", "expressions", "blank")
with gzip.open(os.path.join(EXPR, "features.tsv.gz"), "rt") as f:
    rows = [l.rstrip("\n").split("\t") for l in f]
want = {"ncbi_5333210", "ncbi_113822661", "ncbi_113803873", "ncbi_113808997"}
for r in rows:
    if r[0] in want:
        print("matrix feature:", r)
# which gene_ids exist in GTF and do they have exon rows?
gtf = r"E:\workspace\genome2\genome\GENOME_GTF\GCF_003789085.1.gtf.gz"
targets = {}
for r in rows:
    if r[0] in want:
        targets[r[1]] = r[0]
print("looking for GTF gene ids:", list(targets))
seen = collections.defaultdict(lambda: collections.Counter())
with gzip.open(gtf, "rt") as f:
    for line in f:
        if line[0] == "#": continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9: continue
        m = re.search(r'gene_id "([^"]+)"', p[8])
        if not m: continue
        g = m.group(1)
        if g in targets or ("ncbi_" + g[3:]) in want:
            seen[g][p[2]] += 1
for g, c in seen.items():
    print("GTF", g, dict(c))
