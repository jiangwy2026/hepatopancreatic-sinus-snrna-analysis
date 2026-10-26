import pandas as pd, re, sys
sys.stdout.reconfigure(encoding="utf-8")
P = r"D:\桌面\single cell\文章所需\图5数据\淋巴比对材料.xlsx"
HD = r"D:\桌面\single cell\文章所需\图5数据\human_differ tissue.xlsx"
o = pd.read_excel(P, sheet_name="shrimp-human orthologs"); o.columns = ["human","shrimp"]
h = pd.read_excel(HD); h["g"] = h.gene.astype(str).str.upper().str.strip()
e2g = {}
for _, r in h.iterrows():
    for e in re.findall(r"ENSP\d+", str(r.ENSP_merge)):
        e2g.setdefault(e, set()).add(r.g)
tg = {"JAG1":"113821858","LAMA1":"113802906","THBS4":"113823118","ITGAV":"113821319",
      "VWF":"113824171","KIT":"113804259","MRC1":"113822862"}
for g, gid in tg.items():
    rows = o[o.shrimp.fillna("").astype(str).str.contains(gid, regex=False)]
    for _, r in rows.iterrows():
        es = re.findall(r"ENSP\d+", str(r.human))
        syms = sorted({x for e in es for x in e2g.get(e, [])})
        print(f"{g:6s} gid={gid} row_ensp={len(es)} mapped_symbols={syms[:10]}")
