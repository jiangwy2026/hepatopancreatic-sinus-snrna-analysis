import urllib.request, json, gzip, re, sys
import pandas as pd
sys.stdout.reconfigure(encoding="utf-8")
FASTA = r"D:\桌面\single cell\analysis_output\human_niche_proteins2.fasta"
HFA = r"D:\桌面\single cell\analysis_output\human_refseq_protein.faa.gz"
ORTH = r"E:\workspace\20260729\orthology_project\human_shrimp_all_orthogroups.tsv"
OUT = r"D:\桌面\single cell\analysis_output\step1\sonicparanoid_final_verdict.csv"

genes = [l[1:].strip().split("|") for l in open(FASTA) if l.startswith(">")]
g2np = {}
for g, a in genes:
    try:
        j = json.load(urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{a}.json", timeout=60))
        g2np[g] = sorted({x["id"].split(".")[0] for x in j.get("uniProtKBCrossReferences", [])
                          if x.get("database") == "RefSeq" and x["id"].startswith("NP_")})
    except Exception as e:
        print("FAIL", g, e); g2np[g] = []

present = set()
with gzip.open(HFA, "rt", encoding="utf-8", errors="ignore") as f:
    for l in f:
        if l.startswith(">"):
            present.add(l[1:].split()[0].split(".")[0])

o = pd.read_csv(ORTH, sep="\t", dtype=str).fillna("")
def ids(s):
    return set(x.split(".")[0] for x in re.findall(r"NP_\d+\.\d+", s))
o["hsa"] = o["hsapiens.faa"].map(ids)
o["lva_ids"] = o["pvannamei.faa"].map(lambda s: sorted(set(x.split(".")[0] for x in re.findall(r"XP_\d+\.\d+", s))))

rows = []
for g, nps in g2np.items():
    inp = any(x in present for x in nps)
    hits = o[o["hsa"].map(lambda s: bool(s & set(nps)))]
    if len(hits) == 0:
        rows.append({"human_gene": g, "refseq_np_in_input": inp, "in_orthogroup": False,
                     "group_id": "", "category": "", "shrimp_xp": "", "human_count": 0, "shrimp_count": 0})
    else:
        for _, r in hits.iterrows():
            hc, sc = int(r["human_count"] or 0), int(r["shrimp_count"] or 0)
            cat = "one2one" if hc == 1 and sc == 1 else ("many2many" if hc > 1 and sc > 1 else
                  ("many_hsa_one_shrimp" if hc > 1 else "one_hsa_many_shrimp"))
            rows.append({"human_gene": g, "refseq_np_in_input": inp, "in_orthogroup": True,
                         "group_id": r["group_id"], "category": cat,
                         "shrimp_xp": ";".join(r["lva_ids"]), "human_count": hc, "shrimp_count": sc})
out = pd.DataFrame(rows)
out.to_csv(OUT, index=False, encoding="utf-8-sig")
pd.set_option("display.width", 260); pd.set_option("display.max_colwidth", 60)
print(out.to_string(index=False))
