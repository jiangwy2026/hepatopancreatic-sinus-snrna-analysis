import pandas as pd, re, sys
sys.stdout.reconfigure(encoding="utf-8")
ORTH = r"D:\桌面\single cell\文章所需\图5数据\肝胰腺比对材料.xlsx"
T = r"D:\桌面\single cell\analysis_output\step1\target_ensp.csv"
OUT = r"D:\桌面\single cell\analysis_output\step1\workspace_ortholog_lookup.csv"
orth = pd.read_excel(ORTH, sheet_name="shrimp-human orthologs")
orth["E"] = orth["human"].astype(str).map(lambda s: set(re.findall(r"ENSP\d+", s)))
t = pd.read_csv(T)
rows = []
for gene, grp in t.groupby("gene"):
    es = set(grp["ensp"])
    hits = orth[orth["E"].map(lambda x: bool(x & es))]
    toks = []
    for s in hits["shrimp"]:
        toks.extend([x.strip() for x in str(s).split(",") if x.strip()])
    toks = sorted(set(toks))
    rows.append({"human_gene": gene, "n_ortholog_rows": len(hits),
                 "n_shrimp_tokens": len(toks), "shrimp_tokens": ";".join(toks)})
res = pd.DataFrame(rows).sort_values(["n_ortholog_rows", "human_gene"], ascending=[False, True])
res.to_csv(OUT, index=False, encoding="utf-8-sig")
pd.set_option("display.width", 260); pd.set_option("display.max_colwidth", 120)
print(res.to_string(index=False))
