# -*- coding: utf-8 -*-
import io, sys, time, re
import requests
import xml.etree.ElementTree as ET
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
B = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
PM = ["42000089","42004253","41694335"]
# PubMed -> PMC
print("=== PubMed -> PMC ===")
pmc = {}
for p in PM:
    r = requests.get(B+"elink.fcgi", params={"dbfrom":"pubmed","db":"pmc","id":p,"retmode":"json"}, timeout=30)
    j = r.json(); ids=[]
    for ls in j.get("linksets",[]):
        for db in ls.get("linksetdbs",[]):
            ids += db.get("links",[])
    pmc[p]=ids; print("  PMID %s -> PMC %s" % (p, ids if ids else "(无)")); time.sleep(0.35)
# 完整摘要
print("\n=== 完整摘要 ===")
r = requests.get(B+"efetch.fcgi", params={"db":"pubmed","id":",".join(PM),"retmode":"xml"}, timeout=60)
root = ET.fromstring(r.content)
for art in root.iter("PubmedArticle"):
    pmid = art.findtext(".//PMID")
    title = "".join(art.find(".//ArticleTitle").itertext())
    abst = " ".join("".join(x.itertext()) for x in art.findall(".//AbstractText"))
    print("\n--- PMID %s ---\n%s\n%s" % (pmid, title, abst))
# 若有 PMC，抓全文并找细胞类型定义
for p, ids in pmc.items():
    if not ids: continue
    pid = str(ids[0])
    print("\n=== PMC %s 全文检索（PMID %s）===" % (pid, p))
    try:
        r = requests.get(B+"efetch.fcgi", params={"db":"pmc","id":pid,"retmode":"xml"}, timeout=90)
        txt = re.sub(r"<[^>]+>", " ", r.text)
        txt = re.sub(r"\s+", " ", txt)
        for m in re.finditer(r"([^.]{0,260}(?:B cells|F cells|R cells|E cells|M cells|blister|fibrillar|resorptive)[^.]{0,260}\.)", txt, re.I):
            s = m.group(1).strip()
            if len(s) > 60: print("  •", s[:520])
    except Exception as e:
        print("  全文获取失败:", repr(e)[:120])
