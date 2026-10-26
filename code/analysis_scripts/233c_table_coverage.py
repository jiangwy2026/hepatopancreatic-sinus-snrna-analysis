import pandas as pd, re, sys
sys.stdout.reconfigure(encoding="utf-8")
p = r"D:\桌面\single cell\文章所需\图5数据\肝胰腺比对材料.xlsx"
o = pd.read_excel(p, sheet_name="shrimp-human orthologs")
toks = []
for s in o.shrimp.fillna("").astype(str):
    toks += [x.strip() for x in s.split(",") if x.strip()]
gids = [t.split("|")[0].replace("ncbi_", "") for t in toks]
print("rows", len(o), "shrimp tokens", len(toks), "unique shrimp GeneID", len(set(gids)))
targets = {
 "CXCL12":None,"KITLG":None,"LYVE1":None,"PTN":None,"DAB2":"113808873","KIT":"113804259",
 "CXCR4":"113819525","JAG1":"113821858","DLL1":"113815864","NOTCH1":"113802681",
 "NOTCH2":"113802400","WNT5A":"113806005","FZD5":"113826277","FZD4":"113826275",
 "ANGPT2":"113804527","TEK/TIE1":"113823485","THBS4":"113823118","LAMA1":"113802906",
 "ITGAV":"113821319","ITGB1":"113830322","ITGA4":"113827435","MRC1":"113822862",
 "CSF1R":"113823485","POSTN":"113805211","PROM1":"113802277","VWF":"113824171",
 "CDH5":"113802417","VEGFA":"113804841"}
for g, gid in targets.items():
    hit = [] if gid is None else [t for t, x in zip(toks, gids) if x == gid]
    print(f"{g:9s} {'NA' if gid is None else gid:10s} in_table={len(hit)} {hit[:2]}")
