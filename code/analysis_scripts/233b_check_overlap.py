import pandas as pd, re, sys
sys.stdout.reconfigure(encoding="utf-8")
ORTH = r"D:\桌面\single cell\文章所需\图5数据\肝胰腺比对材料.xlsx"
HD = r"D:\桌面\single cell\文章所需\图5数据\human_differ tissue.xlsx"
o = pd.read_excel(ORTH, sheet_name="shrimp-human orthologs")
h = pd.read_excel(HD)
o["E"] = o["human"].fillna("").astype(str).map(lambda s: set(re.findall(r"ENSP\d+", s)))
h["E"] = h["ENSP_merge"].fillna("").astype(str).map(lambda s: set(re.findall(r"ENSP\d+", s)))
h["g"] = h["gene"].astype(str).str.upper().str.strip()
print("orth rows", len(o), "human genes", h.g.nunique())
for g in ["VEGFA", "MRC1", "DAB2", "ANGPT2", "KIT", "NOTCH1", "JAG1"]:
    if not (h.g == g).any():
        print(f"{g}: not in human_differ"); continue
    es = set().union(*h.loc[h.g == g, "E"].tolist())
    hits = o[o.E.map(lambda x: bool(x & es))]
    toks = [t.strip() for s in hits.shrimp for t in str(s).split(",")]
    print(f"{g}: ensp={len(es)} rows={len(hits)} shrimp={toks[:5]}")
# row-level sample of shared ENSPs
shared = set().union(*o.E.tolist()) & set().union(*h.E.tolist())
print("shared ensp", len(shared))
rows = o[o.E.map(lambda x: bool(x & shared))]
print(rows[["human","shrimp"]].head(3).to_string()[:1500])
