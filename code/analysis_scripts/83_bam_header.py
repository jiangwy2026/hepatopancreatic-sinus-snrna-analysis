# -*- coding: utf-8 -*-
"""Header @PG/@CO + full tag census (all tag names/types) at multiple file positions."""
import zlib, struct, os, collections

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

# --- header
fh = open(BAM, "rb")
gen = bgzf_blocks(fh); buf = bytearray()
def pull(n):
    while len(buf) < n: buf.extend(next(gen))
    out = bytes(buf[:n]); del buf[:n]; return out
pull(4); l_text = struct.unpack("<i", pull(4))[0]; text = pull(l_text).decode(errors="replace")
n_ref = struct.unpack("<i", pull(4))[0]
for i in range(n_ref):
    ln = struct.unpack("<i", pull(4))[0]; pull(ln); pull(4)
fh.close()
print("=== non-@SQ header lines ===")
for line in text.split("\n"):
    if line.startswith("@SQ") or not line.strip(): continue
    print(line[:400])
