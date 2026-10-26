# -*- coding: utf-8 -*-
import io, sys, re, requests
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
B = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
for pid in ("12902777", "13089082"):
    print("\n" + "#"*100); print("PMC", pid); print("#"*100)
    r = requests.get(B+"efetch.fcgi", params={"db":"pmc","id":pid,"retmode":"xml"}, timeout=120)
    txt = re.sub(r"<[^>]+>", " ", r.text); txt = re.sub(r"\s+", " ", txt)
    # 找 Results 里 "cell type / cluster" 相关段落
    for kw in ("cell types", "cell clusters", "cluster identity", "identified the cell", "marker genes",
               "F cell", "R cell", "B cell", "E cell", "M cell"):
        for m in list(re.finditer(re.escape(kw), txt, re.I))[:4]:
            a = max(0, m.start()-300); b = min(len(txt), m.end()+500)
            print("\n--- [%s] ---\n%s" % (kw, txt[a:b].strip()[:800]))
