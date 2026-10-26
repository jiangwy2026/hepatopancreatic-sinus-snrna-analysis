# -*- coding: utf-8 -*-
import io, sys, json, time, re
import requests
import xml.etree.ElementTree as ET
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
B = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
PM = ["42000089","42004253","41694335","41530206","40292183","40246888"]

# 1) PubMed -> GEO
print("=== PubMed -> GEO 关联 ===")
geo_ids = {}
for p in PM:
    try:
        r = requests.get(B + "elink.fcgi", params={"dbfrom":"pubmed","db":"gds","id":p,"retmode":"json"}, timeout=30)
        j = r.json()
        ids = []
        for ls in j.get("linksets", []):
            for db in ls.get("linksetdbs", []):
                if db.get("linkname","").startswith("pubmed_gds"):
                    ids += db.get("links", [])
        geo_ids[p] = ids
        print("  PMID %s -> %s" % (p, ids if ids else "(无)"))
    except Exception as e:
        print("  PMID %s 失败: %r" % (p, e))
    time.sleep(0.35)

# 2) GEO 摘要
allg = sorted({i for v in geo_ids.values() for i in v})
if allg:
    print("\n=== GEO 记录摘要 ===")
    r = requests.get(B + "esummary.fcgi", params={"db":"gds","id":",".join(allg),"retmode":"json"}, timeout=40)
    j = r.json()
    for k, v in j.get("result", {}).items():
        if k == "uids": continue
        print("  GSE%s | %s | %s | n=%s" % (v.get("accession","?"), v.get("gdstype","?"), v.get("title","?")[:150], v.get("n_samples","?")))
        print("     summary:", (v.get("summary","") or "")[:300])

# 3) 摘要
print("\n=== 论文摘要（关键片段）===")
r = requests.get(B + "efetch.fcgi", params={"db":"pubmed","id":",".join(PM),"retmode":"xml"}, timeout=60)
root = ET.fromstring(r.content)
for art in root.iter("PubmedArticle"):
    pmid = art.findtext(".//PMID")
    title = "".join(art.find(".//ArticleTitle").itertext()) if art.find(".//ArticleTitle") is not None else ""
    abst = " ".join("".join(x.itertext()) for x in art.findall(".//AbstractText"))
    print("\n--- PMID %s ---" % pmid)
    print("TITLE:", title)
    print("ABSTRACT:", abst[:1600])
