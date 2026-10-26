import io, sys, requests
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
import xml.etree.ElementTree as ET
B = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
PM = ["42044139","41248710","42580374","42250742","41905082","41638371","41423095","41038398","40865745","40750028","40679295","38323810","37326798","36761772"]
r = requests.get(B+"efetch.fcgi", params={"db":"pubmed","id":",".join(PM),"retmode":"xml"}, timeout=60)
root = ET.fromstring(r.content)
print("=== PubMed 候选（血细胞/淋巴器官单细胞）===")
for art in root.iter("PubmedArticle"):
    pmid = art.findtext(".//PMID")
    t = "".join(art.find(".//ArticleTitle").itertext()) if art.find(".//ArticleTitle") is not None else ""
    jr = art.findtext(".//Journal/ISOAbbreviation") or ""
    yr = art.findtext(".//PubDate/Year") or ""
    ab = " ".join("".join(x.itertext()) for x in art.findall(".//AbstractText"))
    print("\nPMID %s (%s %s)\n  %s" % (pmid, jr, yr, t))
    if ab:
        # 只打印与细胞类型/单细胞相关的片段
        import re
        for m in re.finditer(r"[^.]{0,200}(single-cell|single cell|scRNA|cluster|cell type|subtype)[^.]{0,200}\.", ab, re.I):
            print("   ·", m.group(0).strip()[:300]); break
G = ["200334037","200289112","200273689","200253240","200105787","200061541","200053900","200007456","200004949"]
r = requests.get(B+"esummary.fcgi", params={"db":"gds","id":",".join(G),"retmode":"json"}, timeout=60)
j = r.json()
print("\n\n=== GEO 候选 ===")
for k in G:
    v = j.get("result", {}).get(k)
    if not v: continue
    print("  %s | %s | %s（%s samples）" % (v.get("accession","?"), (v.get("title") or "")[:110], v.get("gdstype","")[:40], v.get("n_samples","?")))
