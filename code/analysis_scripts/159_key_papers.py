import io, sys, re, requests
import xml.etree.ElementTree as ET
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
B = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
for pmid in ("42250742","41638371","42044139","41248710"):
    r = requests.get(B+"efetch.fcgi", params={"db":"pubmed","id":pmid,"retmode":"xml"}, timeout=60)
    root = ET.fromstring(r.content)
    for art in root.iter("PubmedArticle"):
        t = "".join(art.find(".//ArticleTitle").itertext())
        ab = " ".join("".join(x.itertext()) for x in art.findall(".//AbstractText"))
        # 数据可得性/accession
        acc = set(re.findall(r"(GSE\d+|PRJNA\d+|CRA\d+|HRA\d+|OMIX\d+|GSA\d+|SUB\d+)", ab))
        print("="*100); print("PMID", pmid, "|", t); print("="*100)
        print(ab)
        if acc: print("  [数据号]", acc)
        print()
