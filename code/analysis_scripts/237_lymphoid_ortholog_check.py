import pandas as pd, re, sys
sys.stdout.reconfigure(encoding="utf-8")
P = r"D:\桌面\single cell\文章所需\图5数据\淋巴比对材料.xlsx"
o = pd.read_excel(P, sheet_name="shrimp-human orthologs")
o.columns = ["human", "shrimp"]
toks = []
for s in o.shrimp.fillna("").astype(str):
    toks += [x.strip() for x in s.split(",") if x.strip()]
gids = [t.split("|")[0].replace("ncbi_", "") for t in toks]
print("rows", len(o), "tokens", len(toks), "unique gid", len(set(gids)))
tg = {"JAG1":"113821858","LAMA1":"113802906","THBS4":"113823118","ITGAV":"113821319",
      "VWF":"113824171","KIT":"113804259","MRC1":"113822862","DAB2":"113808873",
      "PROM1":"113802277","CDH5":"113802417","POSTN":"113805211","ITGB1":"113830322"}
for g, gid in tg.items():
    hit = [t for t, x in zip(toks, gids) if x == gid]
    print(f"{g:9s} {gid:10s} in_lymphoid_table={len(hit)} {hit[:2]}")
