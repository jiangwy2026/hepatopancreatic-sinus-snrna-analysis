# -*- coding: utf-8 -*-
import io, sys, re, requests
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
B = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
r = requests.get(B+"efetch.fcgi", params={"db":"pmc","id":"13089082","retmode":"xml"}, timeout=120)
txt = re.sub(r"<[^>]+>", " ", r.text); txt = re.sub(r"\s+", " ", txt)
keys = ["trypsin","chymotrypsin","amylase","chitinase","cellulase","cathepsin","carboxypeptidase","lipase",
        "ferritin","apolipoprotein","peritrophin","mucin","vitellogenin","hemocyanin","C-type lectin","integrin"]
print("=== 关键消化/代谢基因在全文中的上下文 ===")
for k in keys:
    ms = list(re.finditer(r"[^.]{0,200}" + re.escape(k) + r"[^.]{0,200}\.", txt, re.I))
    if ms:
        print("\n[%s] %d 处" % (k, len(ms)))
        for m in ms[:2]: print("   •", m.group(0).strip()[:380])
print("\n=== 细胞类型命名的原始描述（前 20 条含 cell type 的句子）===")
seen=set(); n=0
for m in re.finditer(r"[^.]{0,180}(?:B|F|E|M|R) cells?[^.]{0,180}\.", txt):
    s=m.group(0).strip()
    if any(w in s.lower() for w in ("marker","expression","express","specific","cluster","subtype","identified","annotation")):
        if s[:80] in seen: continue
        seen.add(s[:80]); print("  •", s[:330]); n+=1
        if n>=20: break
