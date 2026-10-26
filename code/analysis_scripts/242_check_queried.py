import urllib.request, json, re, sys, pandas as pd
sys.stdout.reconfigure(encoding="utf-8")
T = {"CXCL12":"P48061","KITLG":"P21583","KIT":"P10721","LYVE1":"Q9Y5Y7","PTN":"P21246",
     "CXCR4":"P61073","TEK":"Q02763","TIE1":"P35590","CSF1R":"P07333","CD34":"P28906",
     "NOTCH1":"P46531","JAG1":"P78504","LAMA1":"P25391","THBS4":"P35443","VWF":"P04275"}
o = pd.read_csv(r"E:\workspace\20260729\orthology_project\human_shrimp_all_orthogroups.tsv", sep="\t", dtype=str).fillna("")
human_all = set()
for s in o["hsapiens.faa"]:
    human_all |= set(x.split(".")[0] for x in re.findall(r"NP_\d+\.\d+", s))
print("unique human NP_ in table:", len(human_all))
for g, acc in T.items():
    j = json.load(urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{acc}.json", timeout=60))
    ids = sorted({x["id"].split(".")[0] for x in j.get("uniProtKBCrossReferences", [])
                  if x.get("database") == "RefSeq" and x["id"].startswith("NP_")})
    print(f"{g:8s} uniprot_RefSeq_NP={ids[:5]}  in_table={[x for x in ids if x in human_all]}")
