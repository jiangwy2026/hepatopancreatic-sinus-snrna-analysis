import io, sys, requests
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
B = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
qs = ['Litopenaeus vannamei[All Fields] AND hemocyte AND (single cell OR scRNA)',
      'Litopenaeus vannamei[All Fields] AND (growth immunity trade-off)',
      'Penaeus vannamei[All Fields] AND hemocyte atlas']
for q in qs:
    r = requests.get(B+"esearch.fcgi", params={"db":"sra","term":q,"retmode":"json","retmax":15}, timeout=40)
    ids = r.json().get("esearchresult", {}).get("idlist", [])
    print("=== SRA:", q, "->", len(ids), "条")
    if not ids: continue
    s = requests.get(B+"esummary.fcgi", params={"db":"sra","id":",".join(ids),"retmode":"json"}, timeout=60).json()
    for k, v in s.get("result", {}).items():
        if k == "uids": continue
        print("   %s | %s | %s" % (v.get("accession","?"), v.get("createdate","")[:10], (v.get("title") or "")[:110]))
