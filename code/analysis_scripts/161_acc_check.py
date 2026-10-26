import io, sys, requests
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
B = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
for p in ("42250742","41638371","41248710","42044139","38323810","36761772"):
    r = requests.get(B+"elink.fcgi", params={"dbfrom":"pubmed","db":"gds","id":p,"retmode":"json"}, timeout=30)
    j = r.json(); ids=[]
    for ls in j.get("linksets",[]):
        for db in ls.get("linksetdbs",[]):
            ids += db.get("links",[])
    print("PMID %s -> GEO %s" % (p, ids if ids else "(无)"))
