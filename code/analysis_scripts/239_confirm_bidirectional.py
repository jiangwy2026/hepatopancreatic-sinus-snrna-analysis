import pandas as pd, re, sys
sys.stdout.reconfigure(encoding="utf-8")
P = r"D:\桌面\single cell\文章所需\图5数据\淋巴比对材料.xlsx"
T = r"D:\桌面\single cell\analysis_output\step1\target_ensp.csv"
o = pd.read_excel(P, sheet_name="shrimp-human orthologs"); o.columns = ["human","shrimp"]
t = pd.read_csv(T)
tg = {"JAG1":"113821858","LAMA1":"113802906","THBS4":"113823118","ITGAV":"113821319",
      "VWF":"113824171","KIT":"113804259","MRC1":"113822862","DAB2":"113808873",
      "PROM1":"113802277","CDH5":"113802417","POSTN":"113805211","ITGB1":"113830322"}
for g, gid in tg.items():
    es = set(t.loc[t.gene == g, "ensp"])
    rows = o[o.shrimp.fillna("").astype(str).str.contains(gid, regex=False)]
    hit = False
    for _, r in rows.iterrows():
        rowes = set(re.findall(r"ENSP\d+", str(r.human)))
        common = rowes & es
        if common:
            hit = True
            print(f"{g:6s} gid={gid} CONFIRMED shared_ensp={len(common)} row_ensp={len(rowes)}")
    if not hit:
        print(f"{g:6s} gid={gid} rows={len(rows)} no shared ENSP (target ensp={len(es)})")
