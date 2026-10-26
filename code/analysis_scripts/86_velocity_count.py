# -*- coding: utf-8 -*-
"""
86_velocity_count.py - spliced / unspliced / ambiguous UMI counts from CellRanger 8 BAMs (Windows, no pysam).

Pipeline: zlib-decompress BGZF blocks -> walk BAM records with struct -> classify each read against
GTF exon/intron intervals -> deduplicate molecules by (cell, gene, UMI).

Molecule key (uint64): [2b kind][18b gene][24b umi][20b cell];  kind 0=spliced 1=unspliced 2=ambiguous
Read class: exonic bp == 0 -> unspliced ; exonic bp >= gene-body overlap -> spliced ; else ambiguous
"""
import os, sys, gzip, zlib, struct, time, re, bisect, json
from array import array
import numpy as np
from multiprocessing import Pool

ROOT = r"D:\桌面\single cell"
OUTDIR = os.path.join(ROOT, "analysis_output")
GTF = r"E:\workspace\genome2\genome\GENOME_GTF\GCF_003789085.1.gtf.gz"
EXPR = os.path.join(ROOT, "ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result",
                    "ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result", "1.Expression", "expressions")
BAMDIR = os.path.join(ROOT, "omicsmaster-adrl26080060-on9ab1mqd78xaavc(1)",
                      "10X-3-转录组", "ADRL26080060_order_1")
SAMPLES = [("blank", 0), ("treated", 11753)]
CHUNK = 24 * 1024 * 1024
MARGIN = 256 * 1024
MAXREC = 2_000_000
BSZ = {65: 1, 99: 1, 67: 1, 115: 2, 83: 2, 105: 4, 73: 4, 102: 4}
UMI_TRANS = bytes.maketrans(b"ACGT", b"0123")

def log(*a):
    print(time.strftime("[%H:%M:%S]"), *a, flush=True)

# ================================================================ gene model
def build_model():
    with gzip.open(os.path.join(EXPR, "blank", "features.tsv.gz"), "rt") as f:
        rows = [l.rstrip("\n").split("\t") for l in f]
    feat = [r[0] for r in rows]
    gidx0 = {r[0]: i for i, r in enumerate(rows)}
    gidx1 = {}
    for i, r in enumerate(rows):
        if len(r) > 1 and r[1] and r[1] not in gidx1:
            gidx1[r[1]] = i
    gene_index = {}
    meta = {}; exons = {}
    gid_re = re.compile(r'gene_id "([^"]+)"')
    n_exon = 0
    with gzip.open(GTF, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or (p[2] != "gene" and p[2] != "exon"):
                continue
            m = gid_re.search(p[8])
            if not m:
                continue
            gid = m.group(1)
            gi = gene_index.get(gid)
            if gi is None:
                gi = gidx0.get("ncbi_" + gid[3:]) if gid.startswith("LOC") else None
                if gi is None:
                    gi = gidx1.get(gid)
                if gi is None:
                    gi = gidx0.get(gid)
                if gi is None:
                    continue
                gene_index[gid] = gi
            s = int(p[3]) - 1; e = int(p[4])
            if p[2] == "gene":
                prev = meta.get(gid)
                if prev is None:
                    meta[gid] = [p[0], s, e]
                else:
                    if s < prev[1]: prev[1] = s
                    if e > prev[2]: prev[2] = e
            else:
                exons.setdefault(gid, []).append((s, e)); n_exon += 1
    log("matrix genes %d | GTF genes linked %d | exon rows %d" % (len(feat), len(meta), n_exon))
    contigs = {}; n_noexon = 0
    for gid, (contig, s, e) in meta.items():
        ex = exons.get(gid)
        if not ex:
            n_noexon += 1; continue
        ex.sort()
        merged = []
        for a, b in ex:
            if merged and a <= merged[-1][1]:
                if b > merged[-1][1]: merged[-1][1] = b
            else:
                merged.append([a, b])
        contigs.setdefault(contig, []).append((s, e, gene_index[gid], merged))
    model = {}
    for contig, rows2 in contigs.items():
        rows2.sort()
        n = len(rows2)
        gs = array("i", [0]) * n; ge = array("i", [0]) * n
        gi = array("i", [0]) * n; gm = array("i", [0]) * n
        eoff = array("i", [0]) * (n + 1)
        exs = array("i"); exe = array("i")
        run = 0
        for k, (s, e, idx, merged) in enumerate(rows2):
            gs[k] = s; ge[k] = e; gi[k] = idx
            if e > run: run = e
            gm[k] = run
            for a, b in merged:
                exs.append(a); exe.append(b)
            eoff[k + 1] = len(exs)
        model[contig] = (gs, ge, gi, gm, eoff, exs, exe)
    log("model contigs %d | genes w/o exons skipped %d" % (len(model), n_noexon))
    return feat, model

# ================================================================ BGZF / BAM
def bgzf_blocks(fh):
    while True:
        hdr = fh.read(12)
        if len(hdr) < 12:
            return
        xlen = hdr[10] | (hdr[11] << 8)
        extra = fh.read(xlen)
        if len(extra) < xlen:
            return
        need = (extra[4] | (extra[5] << 8)) + 1 - 12 - xlen - 8
        comp = fh.read(need)
        if len(comp) < need:
            return
        fh.read(8)
        yield zlib.decompress(comp, -15)

def open_bam(path):
    fh = open(path, "rb")
    gen = bgzf_blocks(fh)
    buf = bytearray()
    def pull(n):
        while len(buf) < n:
            buf.extend(next(gen))
        out = bytes(buf[:n]); del buf[:n]; return out
    if pull(4) != b"BAM\x01":
        raise RuntimeError("not a BAM file")
    l_text = struct.unpack("<i", pull(4))[0]
    text = pull(l_text).decode("utf8", "replace")
    n_ref = struct.unpack("<i", pull(4))[0]
    refs = []
    for _ in range(n_ref):
        ln = struct.unpack("<i", pull(4))[0]
        refs.append(pull(ln)[:-1].decode())
        pull(4)
    return fh, gen, bytes(buf), refs, text

def seek_to_frac(fh, frac):
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

def chunk_iter(gen, first_buf, counter, limit=None):
    acc = bytearray(first_buf)
    start = 0; made = 0
    for block in gen:
        acc.extend(block)
        counter[0] += len(block)
        while len(acc) - start >= CHUNK + MARGIN:
            p = start
            while True:
                if p + 4 > len(acc):
                    break
                bs = acc[p] | (acc[p + 1] << 8) | (acc[p + 2] << 16) | (acc[p + 3] << 24)
                if bs > MAXREC or p + 4 + bs > len(acc):
                    break
                p += 4 + bs
                if p - start >= CHUNK:
                    break
            yield bytes(acc[start:p])
            start = p; made += 1
            if limit and made >= limit:
                return
        if start > 2 * CHUNK:
            del acc[:start]; start = 0
    if start < len(acc):
        yield bytes(acc[start:])

# ================================================================ worker
_W = {}
def _init_worker(cell_map, model, refs):
    _W["cell_map"] = cell_map; _W["model"] = model; _W["refs"] = refs

def process_chunk(data):
    cell_map = _W["cell_map"]; model = _W["model"]; refs = _W["refs"]
    keys = []; kap = keys.append
    st = [0] * 11
    L = len(data); p = 0
    up = struct.unpack_from
    cmget = cell_map.get
    nref = len(refs)
    while p + 4 <= L:
        bs = data[p] | (data[p + 1] << 8) | (data[p + 2] << 16) | (data[p + 3] << 24)
        if bs > MAXREC or p + 4 + bs > L:
            break
        o = p + 4; end = o + bs; p = end
        st[0] += 1
        refid, pos, lrn, mapq, bin_, nco, flag, lseq = up("<iiBBHHHi", data, o)
        if flag & 0x904 or nco == 0 or lseq == 0 or mapq < 1:
            st[1] += 1; continue
        contig = refs[refid] if 0 <= refid < nref else None
        if contig is None:
            st[1] += 1; continue
        c = o + 32 + lrn
        t = c + 4 * nco + ((lseq + 1) // 2) + lseq
        cb = None; ur = None
        while t + 3 <= end:
            ty = data[t + 2]
            if ty == 90:
                q = data.index(b"\x00", t + 3, end)
                nm = data[t:t + 2]
                if nm == b"CR":
                    cb = data[t + 3:q]
                    if ur is not None: break
                elif nm == b"UR":
                    ur = data[t + 3:q]
                    if cb is not None: break
                t = q + 1
            elif ty in BSZ:
                t += 3 + BSZ[ty]
            elif ty == 72:
                t = data.index(b"\x00", t + 3, end) + 1
            elif ty == 66:
                t += 8 + up("<i", data, t + 4)[0] * BSZ.get(data[t + 3], 1)
            else:
                break
        if cb is None:
            st[2] += 1; continue
        ci = cmget(cb)
        if ci is None:
            st[2] += 1; continue
        if ur is None or len(ur) != 12:
            st[3] += 1; continue
        gm2 = model.get(contig)
        if gm2 is None:
            st[4] += 1; continue
        blk = []
        if nco == 1:
            op = up("<I", data, c)[0]
            opc = op & 15
            if opc == 0 or opc == 2 or opc == 7 or opc == 8:
                blk.append((pos, pos + (op >> 4)))
        else:
            rp = pos
            for k in range(nco):
                op = up("<I", data, c + 4 * k)[0]
                opc = op & 15
                if opc == 3:
                    rp += op >> 4; continue
                if opc == 1 or opc == 4 or opc == 5 or opc == 6:
                    continue
                blk.append((rp, rp + (op >> 4))); rp += op >> 4
        if not blk:
            st[1] += 1; continue
        gs, ge, gi, gmax, eoff, exs, exe = gm2
        bstart = blk[0][0]; bend = blk[-1][1]
        j = bisect.bisect_left(gs, bend) - 1
        best = -1; best_tot = 0; best_ex = 0; tie = False
        while j >= 0 and gmax[j] > bstart:
            if ge[j] > bstart:
                tot = 0
                for a, b in blk:
                    x = a if a > gs[j] else gs[j]
                    y = b if b < ge[j] else ge[j]
                    if y > x: tot += y - x
                if tot > 0:
                    lo = eoff[j]; hi = eoff[j + 1]; ex = 0
                    for a, b in blk:
                        m2 = bisect.bisect_left(exs, b, lo, hi) - 1
                        while m2 >= lo and exe[m2] > a:
                            xa = exs[m2] if exs[m2] > a else a
                            xb = exe[m2] if exe[m2] < b else b
                            if xb > xa: ex += xb - xa
                            m2 -= 1
                    if tot > best_tot:
                        best_tot = tot; best_ex = ex; best = j; tie = False
                    elif tot == best_tot:
                        tie = True
            j -= 1
        if best < 0:
            st[4] += 1; continue
        if best_ex == 0:
            st[8] += 1
            kind = 1 if not tie else 2
        elif best_ex >= best_tot:
            st[9] += 1
            kind = 0
        else:
            st[10] += 1
            kind = 2
        st[5 + kind] += 1
        try:
            umi = int(ur.translate(UMI_TRANS), 4)
        except ValueError:
            st[3] += 1; continue
        kap((kind << 62) | (gi[best] << 44) | (umi << 20) | ci)
    arr = np.unique(np.array(keys, dtype=np.uint64)) if keys else np.zeros(0, dtype=np.uint64)
    return arr, st

# ================================================================ driver
def run_sample(sample, cell_offset, args):
    with gzip.open(os.path.join(EXPR, sample, "barcodes.tsv.gz"), "rt") as f:
        barcodes = [(x[:-2] if x.endswith("-1") else x) for x in (l.strip() for l in f) if x]
    cell_map = {b.encode(): i + cell_offset for i, b in enumerate(barcodes)}
    bam = os.path.join(BAMDIR, sample, "possorted_genome_bam.bam")
    log("=" * 72)
    log("sample %s | filtered barcodes %d | cell index offset %d" % (sample, len(barcodes), cell_offset))
    feat, model = build_model()
    fh, gen, first_buf, refs, text = open_bam(bam)
    log("bam refs %d | refs in model %d" % (len(refs), sum(1 for r in refs if r in model)))
    if args.get("start_frac"):
        off = seek_to_frac(fh, float(args["start_frac"]))
        log("seeked to %.0f%% of file (%.2f GB)" % (100 * float(args["start_frac"]), off / 2**30))
        gen = bgzf_blocks(fh)
        nb = bytearray()
        while len(nb) < 4 * 1024 * 1024:
            nb.extend(next(gen))
        k = 0
        while k < len(nb) - 8:
            b0 = nb[k] | (nb[k + 1] << 8) | (nb[k + 2] << 16) | (nb[k + 3] << 24)
            if 0 < b0 < MAXREC and k + 4 + b0 <= len(nb) and 0 < nb[k + 12] < 64:
                break
            k += 1
        first_buf = bytes(nb[k:])
        log("resynced after %d bytes" % k)
    NW = int(args.get("workers", 12))
    pool = Pool(NW, initializer=_init_worker, initargs=(cell_map, model, refs))
    counter = [0]; limit = args.get("limit")
    st = [0] * 11; n_chunks = 0
    acc = np.zeros(0, dtype=np.uint64)
    results = []
    t0 = time.time(); t_parent = 0.0; t_wait = 0.0
    def absorb(res):
        nonlocal acc
        arr, s = res
        if len(arr):
            acc = np.concatenate((acc, arr))
        for k2 in range(11): st[k2] += s[k2]
    for chunk in chunk_iter(gen, first_buf, counter, limit=limit):
        ta = time.time()
        results.append(pool.apply_async(process_chunk, (chunk,)))
        n_chunks += 1
        t_parent += time.time() - ta
        while len(results) >= 2 * NW:
            tb = time.time(); absorb(results.pop(0).get()); t_wait += time.time() - tb
        if len(acc) > 60_000_000:
            acc = np.unique(acc)
        if n_chunks % 50 == 0:
            dt = time.time() - t0
            log("  chunks %d | records %.1fM | in-cell %.1fM | decomp %.0f MB/s | %.0fs" % (
                n_chunks, st[0] / 1e6, (st[0] - st[1] - st[2]) / 1e6, counter[0] / 2**20 / dt, dt))
    tb = time.time()
    for r in results:
        absorb(r.get())
    t_wait += time.time() - tb
    pool.close(); pool.join(); fh.close()
    dt = time.time() - t0
    log("parent %.0fs | worker-wait %.0fs" % (t_parent, t_wait))
    acc = np.unique(acc)
    log("done %s: wall %.1fs | records %.1fM | flag_skip %.1fM | no_cell/outside %.1fM | bad_umi %.1fM | no_gene %.1fM" % (
        sample, dt, st[0] / 1e6, st[1] / 1e6, st[2] / 1e6, st[3] / 1e6, st[4] / 1e6))
    log("  read class: fully-intronic %.1fM | fully-exonic %.1fM | exon+intron spanning %.1fM" % (st[8] / 1e6, st[9] / 1e6, st[10] / 1e6))
    log("  kept: spliced %.1fM | unspliced %.1fM | ambiguous %.1fM | molecules %d" % (st[5] / 1e6, st[6] / 1e6, st[7] / 1e6, len(acc)))
    if args.get("dry"):
        return None
    if len(acc) == 0:
        log("no molecules; nothing saved"); return None
    base = acc & np.uint64(0x0FFFFFFFFFFFFFFF)
    kd = (acc >> np.uint64(62)).astype(np.uint8)
    order = np.argsort(base, kind="stable")
    b = base[order]; kd = kd[order]
    same = np.empty(len(b), dtype=bool); same[0] = True
    if len(b) > 1: same[1:] = b[1:] != b[:-1]
    grp = np.cumsum(same)
    multi = np.bincount(grp)[grp] > 1
    n_conf = int(multi.sum())
    kd[multi] = 2
    keep = same | ~multi
    keys = np.unique(b[keep] | (kd[keep].astype(np.uint64) << np.uint64(62)))
    cells = (keys & np.uint64(0xFFFFF)).astype(np.int32)
    genes = ((keys >> np.uint64(44)) & np.uint64(0x3FFFF)).astype(np.int32)
    kinds = (keys >> np.uint64(62)).astype(np.int8)
    out = os.path.join(OUTDIR, "velocity_counts_%s.npz" % sample)
    np.savez_compressed(out, cells=cells, genes=genes, kinds=kinds)
    with open(os.path.join(OUTDIR, "velocity_count_%s.stats.txt" % sample), "w") as f:
        f.write(json.dumps({"sample": sample, "bam": bam, "records": st[0], "flag_skip": st[1],
                            "no_cell_or_outside_filtered": st[2], "bad_umi": st[3], "no_gene": st[4],
                            "spliced_reads": st[5], "unspliced_reads": st[6], "ambiguous_reads": st[7],
                            "fully_intronic": st[8], "fully_exonic": st[9], "spanning": st[10],
                            "molecules_after_conflict": int(len(keys)), "conflict_keys": n_conf,
                            "wall_s": round(dt, 1), "chunks": n_chunks,
                            "decompressed_GB": round(counter[0] / 2**30, 2)}, indent=1))
    log("saved", out)
    return out

if __name__ == "__main__":
    args = {}
    for a in sys.argv[1:]:
        if "=" in a:
            k, v = a.split("=", 1); args[k] = int(v) if v.isdigit() else v
    want = args.get("sample", "blank")
    for s, off in SAMPLES:
        if want in ("all", s):
            run_sample(s, off, args)
