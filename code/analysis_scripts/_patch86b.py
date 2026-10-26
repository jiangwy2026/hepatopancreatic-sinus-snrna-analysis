import re
p = r"D:\桌面\single cell\analysis_output\scripts\86_velocity_count.py"
src = open(p, encoding="utf8").read()

# --- 1. gene id resolution: GTF gene_id "LOC113811791" <-> features col0 "ncbi_113811791"; else symbol match
old = """            g0 = m.group(1)
            gi = gidx1.get(g0)
            if gi is None:
                gi = gidx0.get(g0)
            if gi is None:
                continue
            gid = g0"""
new = """            gid = m.group(1)
            gi = None
            if gid.startswith("LOC"):
                gi = gidx0.get("ncbi_" + gid[3:])
            if gi is None:
                gi = gidx1.get(gid)
            if gi is None:
                gi = gidx0.get(gid)
            if gi is None:
                continue
            gene_index[gid] = gi"""
assert old in src; src = src.replace(old, new)

old = """    meta = {}; exons = {}
    gid_re = re.compile(r'gene_id "([^"]+)"')
    n_exon = 0"""
new = """    meta = {}; exons = {}
    gene_index = {}
    gid_re = re.compile(r'gene_id "([^"]+)"')
    n_exon = 0"""
assert old in src; src = src.replace(old, new)

old = """        contigs.setdefault(contig, []).append((s, e, gidx1.get(gid, gidx0.get(gid)), merged))"""
new = """        contigs.setdefault(contig, []).append((s, e, gene_index[gid], merged))"""
assert old in src; src = src.replace(old, new)

# --- 2. diagnostics: how many reads are fully exonic / spanning / fully intronic
old = """    st = [0, 0, 0, 0, 0, 0, 0, 0]
    L = len(data); p = 0"""
new = """    st = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    L = len(data); p = 0"""
assert old in src; src = src.replace(old, new)

old = """        if best_ex == 0:
            kind = 1 if not tie else 2
        elif best_ex >= best_tot:
            kind = 0
        else:
            kind = 2
        st[5 + kind] += 1"""
new = """        if best_ex == 0:
            st[8] += 1
            kind = 1 if not tie else 2
        elif best_ex >= best_tot:
            st[9] += 1
            kind = 0
        else:
            st[10] += 1
            kind = 2
        st[5 + kind] += 1"""
assert old in src; src = src.replace(old, new)

# --- 3. start_frac: seek into the middle of the BGZF file
old = """def open_bam(path):"""
new = """def seek_to_frac(fh, frac):
    size = os.fstat(fh.fileno()).st_size
    target = int(size * frac); pos = 0
    while pos < target:
        fh.seek(pos)
        hdr = fh.read(12)
        if len(hdr) < 12:
            break
        xlen = hdr[10] | (hdr[11] << 8)
        extra = fh.read(xlen)
        if len(extra) < xlen:
            break
        pos += (extra[4] | (extra[5] << 8)) + 1
    fh.seek(pos)
    return pos

def open_bam(path):"""
assert old in src; src = src.replace(old, new, 1)

old = """    counter = [0]; limit = args.get("limit")"""
new = """    counter = [0]; limit = args.get("limit")
    if args.get("start_frac"):
        off = seek_to_frac(fh, float(args["start_frac"]))
        log("seeked to %.1f%% of file (%.1f GB)" % (100 * float(args["start_frac"]), off / 2**30))
        gen = bgzf_blocks(fh)
        first_buf = b""
        # resync: drop bytes until the first plausible record boundary
        nb = bytearray()
        while len(nb) < 4 * 1024 * 1024:
            nb.extend(next(gen))
        k = 0
        while k < len(nb) - 8:
            b0 = nb[k] | (nb[k+1] << 8) | (nb[k+2] << 16) | (nb[k+3] << 24)
            if 0 < b0 < MAXREC and k + 4 + b0 <= len(nb):
                lrn = nb[k + 12]
                if 0 < lrn < 64:
                    break
            k += 1
        first_buf = bytes(nb[k:])
        log("resynced after", k, "bytes")"""
assert old in src; src = src.replace(old, new)

# --- 4. report diagnostics
old = """    log("  classified reads: spliced %.1fM | unspliced %.1fM | ambiguous %.1fM" % (st[5] / 1e6, st[6] / 1e6, st[7] / 1e6))"""
new = """    log("  classified reads: spliced %.1fM | unspliced %.1fM | ambiguous %.1fM" % (st[5] / 1e6, st[6] / 1e6, st[7] / 1e6))
    log("  overlap diagnostics: fully-exonic %.1fM | exon+intron spanning %.1fM | fully-intronic %.1fM" % (st[9] / 1e6, st[10] / 1e6, st[8] / 1e6))"""
assert old in src; src = src.replace(old, new)

old = """                            "chunks": n_chunks,"""
new = """                            "chunks": n_chunks, "fully_exonic": st[9], "spanning": st[10], "fully_intronic": st[8],"""
assert old in src; src = src.replace(old, new)

open(p, "w", encoding="utf8").write(src)
print("patched ok")
