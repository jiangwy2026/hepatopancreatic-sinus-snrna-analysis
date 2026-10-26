import io, sys, json, requests
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
# 1) Europe PMC
r = requests.get("https://www.ebi.ac.uk/europepmc/webservices/rest/search",
                 params={"query": "EXT_ID:42250742 AND SRC:MED", "format": "json", "resultType": "core"}, timeout=40)
j = r.json()
res = j.get("resultList", {}).get("result", [])
if res:
    x = res[0]
    print("=== Europe PMC ===")
    for k in ("title","journalInfo","doi","pmcid","isOpenAccess","inEPMC","hasSuppl","hasTextMinedTerms","hasPDF"):
        print("  %-18s %s" % (k, str(x.get(k))[:200]))
    print("  abstract head:", (x.get("abstractText") or "")[:200])
else:
    print("Europe PMC 无记录")
print()
# 2) Crossref DOI
r = requests.get("https://api.crossref.org/works", params={"query.bibliographic": "Single-cell transcriptomics reveals the hemocyte atlas and molecular mechanisms underlying the growth-immunity trade-off in Litopenaeus vannamei", "rows": 2}, timeout=40)
j = r.json()
for it in j.get("message", {}).get("items", [])[:2]:
    print("=== Crossref ===")
    print("  DOI:", it.get("DOI"), "| journal:", it.get("container-title"), "| vol/page:", it.get("volume"), it.get("page"))
