# -*- coding: utf-8 -*-
"""Diagnostic: dump 5 raw records from mid-file, then census CB presence at several file positions."""
import zlib, struct, os, sys, time, collections

BAM = r"D:\桌面\single cell\omicsmaster-adrl26080060-on9ab1mqd78xaavc(1)\10X-3-转录组\ADRL26080060_order_1\blank\possorted_genome_bam.bam"
SIZE = os.path.getsize(BAM)

def bgzf_blocks(fh):
    while True:
        hdr = fh.read(12)
        if len(hdr) < 12: return
        xlen = int.from_bytes(hdr[10:12], "little")
        extra = fh.read(xlen)
        bsize = int.from_bytes(extra[4:6], "little") + 1
        comp = fh.read(bsize - 12 - xlen - 8); fh.read(8)
        yield zlib.decompress(comp, -15)

def seek_block(fh, frac):
    target = int(SIZE * frac); pos = 0
    while pos < target:
        fh.seek(pos); hdr = fh.read(12)
        if len(hdr) < 12: break
        xlen = int.from_bytes(hdr[10:12], "little")
        extra = fh.read(xlen)
        if len(extra) < xlen: break
        pos += int.from_bytes(extra[4:6], "little") + 1
    fh.seek(pos); return pos

def records_from(fh, nblocks, resync_drop=2):
    gen = bgzf_blocks(fh)
    carry = b""; out = []; nb = 0
    first = True
    for d in gen:
        nb += 1
        data = carry + d
        L = len(data); p = 0
        while p + 4 <= L:
            bs = int.from_bytes(data[p:p+4], "little")
            if bs > 200000: break
            if p + 4 + bs > L: break
            out.append(bytes(data[p+4:p+4+bs])); p += 4 + bs
        carry = data[p:]
        if nb >= nblocks: break
    return out[nblocks:] if False else out

fh = open(BAM, "rb")
pos = seek_block(fh, 0.4)
recs = records_from(fh, 3)
print("=== 5 raw records at 40%% (offset %.1f GB, block %d) ===" % (pos/2**30, 0))
for r in recs[3:8]:
    lrn, mapq, bin_, nco, flag, lseq = struct.unpack_from("<BBHHHi", r, 8)
    p = 32 + lrn
    cig = struct.unpack_from("<%dI" % nco, r, p) if nco else ()
    cigs = "".join("%d%s" % (c >> 4, "MIDNSHP=X"[c & 15]) for c in cig)
    name = r[32:32+lrn-1].decode(errors="replace")
    tagstart = p + 4*nco + ((lseq+1)//2) + lseq
    print("flag=%d lseq=%d nco=%d namelen=%d recsize=%d name=%s" % (flag, lseq, nco, lrn, len(r)+4, name))
    print("   cigar:", cigs[:80])
    print("   tag bytes:", r[tagstart:tagstart+90].hex())
    print("   tag ascii:", repr(r[tagstart:tagstart+90]))
fh.close()

print()
print("=== CB census across file positions ===")
_TAGFIX = {"A":1,"c":1,"C":1,"s":2,"S":2,"i":4,"I":4,"f":4}
def has_tags(r):
    lrn, mapq, bin_, nco, flag, lseq = struct.unpack_from("<BBHHHi", r, 8)
    p = 32 + lrn + 4 * nco + ((lseq + 1) // 2) + lseq
    found = {}
    end = len(r)
    while p + 3 <= end:
        t = r[p:p+2]; ty = r[p+2]
        if ty in _TAGFIX: p += 3 + _TAGFIX[ty]
        elif ty in (90, 72):
            q = r.index(b"\x00", p + 3)
            if t in (b"CB", b"UB", b"RE", b"GX"): found[t.decode()] = r[p+3:q].decode(errors="replace")
            p = q + 1
        elif ty == 66:
            sub = r[p+3]; n = struct.unpack_from("<i", r, p+4)[0]; p += 8 + n*_TAGFIX[chr(sub)]
        else: break
    return found

for frac in (0.02, 0.10, 0.25, 0.40, 0.55, 0.70, 0.85, 0.97):
    fh = open(BAM, "rb"); seek_block(fh, frac)
    recs = records_from(fh, 6)
    fh.close()
    c = collections.Counter(); ex = None
    for r in recs:
        try: tg = has_tags(r)
        except Exception: continue
        if tg: c["tagged"] += 1
        if "CB" in tg:
            c["CB"] += 1
            if ex is None: ex = tg
        if "RE" in tg: c["RE"] += 1
        if "GX" in tg: c["GX"] += 1
    print("frac %.2f  n=%d  tagged=%d  CB=%d  RE=%d  GX=%d  example=%s" % (frac, len(recs), c["tagged"], c["CB"], c["RE"], c["GX"], ex))
