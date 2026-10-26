# -*- coding: utf-8 -*-
"""Tag census v2 with a chunk-based parser (bytes + absolute offsets): what tags do in-cell reads carry?"""
import zlib, struct, os, collections, time

BAM = r"D:\桌面\single cell\omicsmaster-adrl26080060-on9ab1mqd78xaavc(1)\10X-3-转录组\ADRL26080060_order_1\blank\possorted_genome_bam.bam"
SIZE = os.path.getsize(BAM)
U16 = struct.Struct("<BBHHHi").unpack_from

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

_TF = {65:1, 99:1, 67:1, 115:2, 83:2, 105:4, 73:4, 102:4}
def walk_chunk(data, base, cb_count, re_dist, gx_cnt, taginv, nco_hist, example, n_lim=None):
    L = len(data); p = 0; n = 0
    while p + 4 <= L:
        bs = data[p] | (data[p+1] << 8) | (data[p+2] << 16) | (data[p+3] << 24)
        if bs > 500000 or p + 4 + bs > L: break
        o = p + 4
        lrn, mapq, bin_, nco, flag, lseq = U16(data, o + 8)
        nco_hist[min(nco, 6)] += 1
        t = o + 32 + lrn + 4 * nco + ((lseq + 1) // 2) + lseq
        end = o + bs
        cb = None; tag_here = []
        while t + 3 <= end:
            ty = chr(data[t+2]); name = data[t:t+2]
            if ty in "AcCsSiIf":
                ln = _TF[ord(ty)]
                if name in (b"RE", b"GN"):
                    tag_here.append((name.decode(), data[t+3:t+3+ln].decode(errors="replace")))
                t += 3 + ln
            elif ty in "ZH":
                q = data.index(b"\x00", t + 3, end)
                v = data[t+3:q]
                if name == b"CR": cb = v
                elif name in (b"GX", b"GN", b"RE", b"UR"):
                    tag_here.append((name.decode(), v.decode(errors="replace")))
                t = q + 1
            elif ty == "B":
                sub = data[t+3]; cnt = struct.unpack_from("<i", data, t+4)[0]
                t += 8 + cnt * _TF.get(sub, 1)
            else: break
        # tag inventory (first 200 records only)
        if n < 200:
            tt = o + 32 + lrn + 4 * nco + ((lseq + 1) // 2) + lseq
            while tt + 3 <= end:
                ty2 = chr(data[tt+2]); nm2 = data[tt:tt+2].decode(errors="replace")
                taginv["%s:%s" % (nm2, ty2)] += 1
                if ty2 in "AcCsSiIf": tt += 3 + _TF[ord(ty2)]
                elif ty2 in "ZH":
                    q = data.index(b"\x00", tt + 3, end); tt = q + 1
                elif ty2 == "B":
                    sub = data[tt+3]; cnt = struct.unpack_from("<i", data, tt+4)[0]; tt += 8 + cnt * _TF.get(sub, 1)
                else: break
        if cb is not None:
            cb_count[0] += 1
            re_v = dict(tag_here).get("RE")
            if re_v is not None: re_dist[re_v] += 1
            if "GX" in dict(tag_here): gx_cnt[0] += 1
            if example[0] is None and "GX" in dict(tag_here):
                example[0] = tag_here
        p += 4 + bs; n += 1
    return p, n

taginv = collections.Counter(); re_dist = collections.Counter(); nco_hist = collections.Counter()
cb_count = [0]; gx_cnt = [0]; example = [None]
tot = 0; t0 = time.time()
for frac in (0.05, 0.2, 0.35, 0.5, 0.65, 0.8, 0.95):
    fh = open(BAM, "rb"); seek_block(fh, frac)
    gen = bgzf_blocks(fh); carry = b""; nb = 0
    while nb < 20:
        try: d = next(gen)
        except StopIteration: break
        nb += 1
        data = carry + d
        p, n = walk_chunk(data, 0, cb_count, re_dist, gx_cnt, taginv, nco_hist, example)
        carry = data[p:]; tot += n
    fh.close()
dt = time.time() - t0
print("scanned %d records in %.1fs -> %.0f rec/s (parser incl. decompress)" % (tot, dt, tot/dt))
print("records with CR tag: %d (%.1f%%)" % (cb_count[0], 100*cb_count[0]/tot))
print("records with CR+GX   : %d (%.1f%% of total, %.1f%% of CR)" % (gx_cnt[0], 100*gx_cnt[0]/tot, 100*gx_cnt[0]/max(cb_count[0],1)))
print("RE dist among CR reads:", re_dist.most_common())
print("n_cigar_op histogram (capped 6+):", nco_hist.most_common())
print("tag:type inventory (first 200 recs/block):", taginv.most_common(20))
print("example in-cell tags:", example[0])
