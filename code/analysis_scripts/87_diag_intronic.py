# -*- coding: utf-8 -*-
"""87_diag_intronic.py - why are there no unspliced reads? Inspect RE-tagged non-exonic reads."""
import importlib.util, os, gzip, bisect, struct, collections, time
spec = importlib.util.spec_from_file_location("vc", r"D:\桌面\single cell\analysis_output\scripts\86_velocity_count.py")
vc = importlib.util.module_from_spec(spec); spec.loader.exec_module(vc)

feat, model = vc.build_model()
sample = "blank"
with gzip.open(os.path.join(vc.EXPR, sample, "barcodes.tsv.gz"), "rt") as f:
    barcodes = [(x[:-2] if x.endswith("-1") else x) for x in (l.strip() for l in f) if x]
cell_map = {b.encode(): i for i, b in enumerate(barcodes)}
bam = os.path.join(vc.BAMDIR, sample, "possorted_genome_bam.bam")
fh, gen, first_buf, refs, text = vc.open_bam(bam)
off = vc.seek_to_frac(fh, 0.7)
gen = vc.bgzf_blocks(fh)
nb = bytearray()
while len(nb) < 4 * 1024 * 1024:
    nb.extend(next(gen))
first_buf = bytes(nb)
data = first_buf

up = struct.unpack_from
BSZ = vc.BSZ
buckets = collections.Counter()
re_by_class = collections.Counter()
examples = collections.defaultdict(list)
contig_hits = collections.Counter()
n = 0; t0 = time.time()
L = len(data); p = 0
while p + 4 <= L:
    bs = data[p] | (data[p+1] << 8) | (data[p+2] << 16) | (data[p+3] << 24)
    if bs > vc.MAXREC or p + 4 + bs > L: break
    o = p + 4; end = o + bs; p = end
    refid, pos, lrn, mapq, bin_, nco, flag, lseq = up("<iiBBHHHi", data, o)
    if flag & 0x904 or nco == 0 or lseq == 0 or mapq < 1: continue
    contig = refs[refid]
    c = o + 32 + lrn
    t = c + 4*nco + ((lseq+1)//2) + lseq
    cb = None; ur = None; re_ = None; gx = None
    while t + 3 <= end:
        ty = data[t+2]
        if ty == 90:
            q = data.index(b"\x00", t+3, end); nm = data[t:t+2]
            if nm == b"CR": cb = data[t+3:q]
            elif nm == b"UR": ur = data[t+3:q]
            elif nm == b"RE": re_ = data[t+3:q]
            elif nm == b"GX": gx = data[t+3:q]
            t = q + 1
        elif ty in BSZ: t += 3 + BSZ[ty]
        elif ty == 72: t = data.index(b"\x00", t+3, end) + 1
        elif ty == 66: t += 8 + up("<i", data, t+4)[0] * BSZ.get(data[t+3], 1)
        else: break
    if cb is None or cb not in cell_map: continue
    n += 1
    contig_hits[contig] += 1
    blk = []
    rp = pos
    for k in range(nco):
        op = up("<I", data, c + 4*k)[0]; opc = op & 15
        if opc == 3: rp += op >> 4; continue
        if opc in (1, 4, 5, 6): continue
        blk.append((rp, rp + (op >> 4))); rp += op >> 4
    if not blk: continue
    gm2 = model.get(contig)
    rex = re_.decode() if re_ else "?"
    if gm2 is None:
        buckets[("no_model_contig", rex)] += 1
        continue
    gs, ge, gi, gmax, eoff, exs, exe = gm2
    bstart = blk[0][0]; bend = blk[-1][1]
    j = bisect.bisect_left(gs, bend) - 1
    best = -1; best_tot = 0; best_ex = 0
    while j >= 0 and gmax[j] > bstart:
        if ge[j] > bstart:
            tot = 0
            for a, b in blk:
                x = a if a > gs[j] else gs[j]; y = b if b < ge[j] else ge[j]
                if y > x: tot += y - x
            if tot > 0:
                lo = eoff[j]; hi = eoff[j+1]; ex = 0
                for a, b in blk:
                    m2 = bisect.bisect_left(exs, b, lo, hi) - 1
                    while m2 >= lo and exe[m2] > a:
                        xa = exs[m2] if exs[m2] > a else a
                        xb = exe[m2] if exe[m2] < b else b
                        if xb > xa: ex += xb - xa
                        m2 -= 1
                if tot > best_tot:
                    best_tot = tot; best_ex = ex; best = j
        j -= 1
    alen = sum(b - a for a, b in blk)
    if best < 0:
        cls = "no_gene_overlap"
    elif best_ex == 0:
        cls = "fully_intronic"
    elif best_ex >= best_tot:
        cls = "fully_exonic"
    else:
        cls = "spanning"
    buckets[(cls, rex)] += 1
    re_by_class[(cls, bool(gx))] += 1
    if cls != "fully_exonic" and len(examples[cls]) < 6:
        examples[cls].append((contig, pos, blk, best_tot, best_ex, alen, rex,
                             gx.decode() if gx else None))
dt = time.time() - t0
print("in-cell reads inspected: %d  (%.1fs, %.0f rec/s)" % (n, dt, n/dt))
print("top contigs:", contig_hits.most_common(5))
print("--- class x RE tag ---")
for k, v in buckets.most_common(20):
    print("   %-18s RE=%-3s %8d (%.1f%%)" % (k[0], k[1], v, 100*v/n))
print("--- class x has_GX ---")
for k, v in re_by_class.most_common(20):
    print("   %-18s GX=%-5s %8d (%.1f%%)" % (k[0], k[1], v, 100*v/n))
print("--- examples of non-fully-exonic reads ---")
for cls, exs_ in examples.items():
    for e in exs_[:3]:
        print("  %s | contig %s pos %d blocks %s gene_overlap %d exon_overlap %d aligned %d RE=%s GX=%s" % (
            cls, e[0], e[1], e[2][:3], e[3], e[4], e[5], e[6], e[7]))
