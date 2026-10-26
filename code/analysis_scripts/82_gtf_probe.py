import gzip, collections, re, time
gtf = r"E:\workspace\genome2\genome\GENOME_GTF\GCF_003789085.1.gtf.gz"
import os
print("gtf exists:", os.path.exists(gtf), "size MB:", round(os.path.getsize(gtf)/2**20,1) if os.path.exists(gtf) else None)
t0=time.time()
ftypes = collections.Counter(); chroms=collections.Counter(); genes=set(); chroms6=0
exon_rows=0
with gzip.open(gtf, "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("#"): continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9: continue
        c, src, ft, st, en, sc, st_, fr, attrs = p
        ftypes[ft]+=1; chroms[c]+=1
        if ft == "exon": exon_rows += 1
        if ft == "gene":
            m = re.search(r'gene_id "([^"]+)"', attrs)
            if m: genes.add(m.group(1))
print("elapsed", round(time.time()-t0,1), "s")
print("feature types:", ftypes.most_common(15))
print("n chroms:", len(chroms), "top:", chroms.most_common(5))
print("n genes:", len(genes))
print("exon rows:", exon_rows)
