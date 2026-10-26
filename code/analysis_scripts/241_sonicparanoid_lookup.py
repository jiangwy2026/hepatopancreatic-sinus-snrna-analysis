import urllib.request, json, re, sys, gzip
import pandas as pd
sys.stdout.reconfigure(encoding="utf-8")

FASTA = r"D:\桌面\single cell\analysis_output\human_niche_proteins2.fasta"
ORTH = r"E:\workspace\20260729\orthology_project\human_shrimp_all_orthogroups.tsv"
FEAT = r"D:\桌面\single cell\analysis_output\GCF_003789085.1_feature_table.txt.gz"
OUT = r"D:\桌面\single cell\analysis_output\step1\sonicparanoid_niche_lookup.csv"

genes = [l[1:].strip().split("|") for l in open(FASTA) if l.startswith(">")]
# 1) human RefSeq NP_ accessions per gene
g2np = {}
for g, a in genes:
    try:
        j = json.load(urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{a}.json", timeout=60))
        ids = [x["id"].split(".")[0] for x in j.get("uniProtKBCrossReferences", [])
               if x.get("database") == "RefSeq" and x["id"].startswith("NP_")]
        g2np[g] = set(ids)
    except Exception as e:
        print("FAIL", g, e); g2np[g] = set()
print("genes with NP_ ids:", sum(1 for v in g2np.values() if v))

# 2) orthogroup table
o = pd.read_csv(ORTH, sep="\t", dtype=str).fillna("")
def ids(s, pat):
    return set(x.split(".")[0] for x in re.findall(pat, s))
o["hsa"] = o["hsapiens.faa"].map(lambda s: ids(s, r"NP_\d+\.\d+"))
o["lva"] = o["pvannamei.faa"].map(lambda s: ids(s, r"XP_\d+\.\d+"))

# 3) shrimp XP_ -> GeneID
rows = []
with gzip.open(FEAT, "rt", encoding="utf-8", errors="ignore") as f:
    hdr = f.readline().lstrip("#").strip().split("\t")
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < len(hdr): continue
        d = dict(zip(hdr, p))
        if d.get("feature") != "mRNA": continue
        rows.append((d.get("related_accession", "").split(".")[0], d.get("GeneID", "")))
xp2gid = {}
for xp, gid in rows:
    if xp: xp2gid.setdefault(xp, gid)

res = []
for g, nps in g2np.items():
    if not nps:
        res.append({"human_gene": g, "found": False, "group_id": "", "category": "",
                    "shrimp_gene_ids": "", "shrimp_xp": "", "n_hsa_ids": 0, "n_lva_ids": 0})
        continue
    hits = o[o["hsa"].map(lambda s: bool(s & nps))]
    if len(hits) == 0:
        res.append({"human_gene": g, "found": False, "group_id": "", "category": "",
                    "shrimp_gene_ids": "", "shrimp_xp": "", "n_hsa_ids": 0, "n_lva_ids": 0})
        continue
    for _, r in hits.iterrows():
        lva = sorted(r["lva"])
        gids = sorted({xp2gid.get(x, "") for x in lva} - {""})
        hc, sc = int(r["human_count"] or 0), int(r["shrimp_count"] or 0)
        cat = "one2one" if hc == 1 and sc == 1 else ("many2many" if hc > 1 and sc > 1 else
              ("many_hsa_one_shrimp" if hc > 1 else "one_hsa_many_shrimp"))
        res.append({"human_gene": g, "found": True, "group_id": r["group_id"], "category": cat,
                    "shrimp_gene_ids": ";".join(gids), "shrimp_xp": ";".join(lva),
                    "n_hsa_ids": hc, "n_lva_ids": sc, "shared_np": ";".join(sorted(nps & r["hsa"]))})
out = pd.DataFrame(res)
out.to_csv(OUT, index=False, encoding="utf-8-sig")
pd.set_option("display.width", 300); pd.set_option("display.max_colwidth", 90)
print(out.to_string(index=False))
