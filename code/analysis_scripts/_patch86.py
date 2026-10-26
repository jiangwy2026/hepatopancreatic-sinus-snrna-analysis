import re
p = r"D:\桌面\single cell\analysis_output\scripts\86_velocity_count.py"
src = open(p, encoding="utf8").read()

old = """    with gzip.open(os.path.join(EXPR, "blank", "features.tsv.gz"), "rt") as f:
        feat = [l.rstrip("\\n").split("\\t")[0] for l in f]
    gidx = {g: i for i, g in enumerate(feat)}
    meta = {}; exons = {}
    gid_re = re.compile(r'gene_id "([^"]+)"')
    n_exon = 0"""
new = """    with gzip.open(os.path.join(EXPR, "blank", "features.tsv.gz"), "rt") as f:
        rows = [l.rstrip("\\n").split("\\t") for l in f]
    feat = [r[0] for r in rows]
    gidx0 = {r[0]: i for i, r in enumerate(rows)}
    gidx1 = {r[1]: i for i, r in enumerate(rows) if len(r) > 1 and r[1] not in gidx0 or True}
    gidx1 = {}
    for i, r in enumerate(rows):
        if len(r) > 1 and r[1] and r[1] not in gidx1:
            gidx1[r[1]] = i
    gidx = gidx0
    meta = {}; exons = {}
    gid_re = re.compile(r'gene_id "([^"]+)"')
    n_exon = 0"""
assert old in src
src = src.replace(old, new)

old2 = """            gid = "ncbi_" + m.group(1)
            if gid not in gidx:
                continue"""
new2 = """            g0 = m.group(1)
            gi = gidx1.get(g0)
            if gi is None:
                gi = gidx0.get(g0)
            if gi is None:
                continue
            gid = g0"""
assert old2 in src
src = src.replace(old2, new2)

old3 = """        contigs.setdefault(contig, []).append((s, e, gidx[gid], merged))"""
new3 = """        contigs.setdefault(contig, []).append((s, e, gidx1.get(gid, gidx0.get(gid)), merged))"""
assert old3 in src
src = src.replace(old3, new3)

old4 = """    feat = [r[0] for r in rows]"""
# (already replaced above)

# guard empty key array
old5 = """    base = acc & np.uint64(0x0FFFFFFFFFFFFFFF)"""
new5 = """    if len(acc) == 0:
        log("no molecules - aborting save"); return None
    base = acc & np.uint64(0x0FFFFFFFFFFFFFFF)"""
assert old5 in src
src = src.replace(old5, new5)

# instrumentation: time in parent (decompress+walk+dispatch) vs waiting
old6 = """    for chunk in chunk_iter(gen, first_buf, counter, limit=limit):
        results.append(pool.apply_async(process_chunk, (chunk,)))
        n_chunks += 1
        while len(results) >= 2 * NW:
            absorb(results.pop(0).get())"""
new6 = """    t_parent = 0.0; t_wait = 0.0
    for chunk in chunk_iter(gen, first_buf, counter, limit=limit):
        ta = time.time()
        results.append(pool.apply_async(process_chunk, (chunk,)))
        n_chunks += 1
        t_parent += time.time() - ta
        while len(results) >= 2 * NW:
            tb = time.time()
            absorb(results.pop(0).get())
            t_wait += time.time() - tb"""
assert old6 in src
src = src.replace(old6, new6)

old7 = """    for r in results:
        absorb(r.get())
    pool.close(); pool.join(); fh.close()"""
new7 = """    tb = time.time()
    for r in results:
        absorb(r.get())
    t_wait += time.time() - tb
    pool.close(); pool.join(); fh.close()
    log("  timing: parent(dispatch) %.0fs | waiting for workers %.0fs" % (t_parent, t_wait))"""
assert old7 in src
src = src.replace(old7, new7)

open(p, "w", encoding="utf8").write(src)
print("patched")
