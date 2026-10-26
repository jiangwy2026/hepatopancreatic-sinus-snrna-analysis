# -*- coding: utf-8 -*-
"""157_lookup_C0_genes.py - 用 NCBI Gene 查 C0 特异基因的描述"""
import io, sys, time, json, requests
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
B = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
ids = ["113819543","113823493","113830563","113830564","113830566","113815249","113819837","113827561","113805372","113810905","113813598","113806311","113817924"]
r = requests.get(B+"esummary.fcgi", params={"db":"gene","id":",".join(ids),"retmode":"json"}, timeout=60)
j = r.json()
res = j.get("result", {})
print("=== NCBI Gene 注释 ===")
for gid in ids:
    v = res.get(gid)
    if not v:
        print("  LOC%s : (无记录)" % gid); continue
    print("  LOC%-12s | %s | %s | %s" % (gid, v.get("name","?"), v.get("description","?"), v.get("summary","")[:220]))
