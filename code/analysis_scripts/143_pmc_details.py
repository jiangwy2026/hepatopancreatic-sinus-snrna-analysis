# -*- coding: utf-8 -*-
import io, sys, re, requests
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
B = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
for pid, tag in (("13089082", "Animal Nutrition 2026 (11 cell types)"), ("12902777", "Front Immunol 2026 (GSE306674)")):
    print("\n" + "="*100); print("PMC", pid, "|", tag); print("="*100)
    r = requests.get(B+"efetch.fcgi", params={"db":"pmc","id":pid,"retmode":"xml"}, timeout=120)
    txt = re.sub(r"<[^>]+>", " ", r.text); txt = re.sub(r"\s+", " ", txt)
    # 数据可得性
    for pat in (r"Data availability[^.]{0,600}\.", r"(GSE\d+|PRJNA\d+|CRA\d+|CNP\d+|HRA\d+|OMIX\d+)", r"accession[^.]{0,300}\."):
        for m in list(re.finditer(pat, txt, re.I))[:6]:
            s = m.group(0).strip()
            if len(s) > 20: print("  [数据]", s[:400]); break
    # 各细胞类型 marker 描述
    for ct in ("B cell", "F cell", "E cell", "R cell", "M cell", "B cells", "F cells", "E cells", "R cells", "M cells"):
        for m in list(re.finditer(r"[^.]{0,220}" + ct + r"[^.]{0,220}\.", txt, re.I))[:3]:
            s = m.group(0).strip()
            if any(k in s.lower() for k in ("marker", "express", "enrich", "character", "identif", "defined", "gene")):
                print("  [%s]" % ct, s[:420])
