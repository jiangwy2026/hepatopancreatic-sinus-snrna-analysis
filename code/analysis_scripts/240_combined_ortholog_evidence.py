import pandas as pd, re, sys
sys.stdout.reconfigure(encoding="utf-8")
HP = r"D:\桌面\single cell\文章所需\图5数据\肝胰腺比对材料.xlsx"
LY = r"D:\桌面\single cell\文章所需\图5数据\淋巴比对材料.xlsx"
BL = r"D:\桌面\single cell\analysis_output\step1\human_shrimp_niche_orthologs_full.csv"
hp_in = pd.read_excel(HP, sheet_name="tpm>100")
ly_in = pd.read_excel(LY, sheet_name="Sheet1")
hp_o = pd.read_excel(HP, sheet_name="shrimp-human orthologs"); hp_o.columns=["human","shrimp"]
ly_o = pd.read_excel(LY, sheet_name="shrimp-human orthologs"); ly_o.columns=["human","shrimp"]
def gids(df):
    out=set()
    for s in df.shrimp.fillna("").astype(str):
        for t in s.split(","):
            t=t.strip()
            if t: out.add(t.split("|")[0].replace("ncbi_",""))
    return out
HPG, LYG = gids(hp_o), gids(ly_o)
# search input gene lists by symbol/description
def find_in(df, cols, pats):
    txt = df[cols].fillna("").astype(str).agg(" | ".join, axis=1).str.lower()
    for p in pats:
        m = txt.str.contains(p, regex=False)
        if m.any():
            print(f"   input match '{p}':", df.loc[m, cols[:2]].head(3).values.tolist())
print("== hepatopancreas input (983 genes) search ==")
find_in(hp_in, ["Symbol","Description"], ["cxcl","kitlg","stem cell factor","lyve","pleiotrophin","pty"])
print("== lymphoid input (800 genes) search ==")
find_in(ly_in, ["Gene Name","Gene ID"], ["cxcl","kitlg","lyve","pleiotrophin","ptn"])
# blast table
b = pd.read_csv(BL)
bm = dict(zip(b.human_gene, b.verdict))
def bl(g):
    return bm.get(g, "no hit")
rows = []
targets = ["CXCL12","KITLG","LYVE1","PTN","DAB2","KIT","CXCR4","JAG1","DLL1","NOTCH1","NOTCH2",
           "WNT5A","FZD5","FZD4","ANGPT2","TEK","TIE1","THBS4","LAMA1","ITGAV","ITGB1","ITGA4",
           "MRC1","CSF1R","POSTN","PROM1","VWF","CDH5","VEGFA"]
gmap = {"DAB2":"113808873","KIT":"113804259","CXCR4":"113819525","JAG1":"113821858","DLL1":"113815864",
        "NOTCH1":"113802681","NOTCH2":"113802400","WNT5A":"113806005","FZD5":"113826277","FZD4":"113826275",
        "ANGPT2":"113804527","TEK":"113823485","TIE1":"113823485","THBS4":"113823118","LAMA1":"113802906",
        "ITGAV":"113821319","ITGB1":"113830322","ITGA4":"113827435","MRC1":"113822862","CSF1R":"113823485",
        "POSTN":"113805211","PROM1":"113802277","VWF":"113824171","CDH5":"113802417","VEGFA":"113804841"}
for g in targets:
    gid = gmap.get(g)
    rows.append({"human_gene":g,"shrimp_gene_id":gid or "",
                 "in_hepatopancreas_table": (gid in HPG) if gid else False,
                 "in_lymphoid_table": (gid in LYG) if gid else False,
                 "blastp_verdict": bl(g)})
out = pd.DataFrame(rows)
out.to_csv(r"D:\桌面\single cell\analysis_output\step1\combined_ortholog_evidence.csv", index=False, encoding="utf-8-sig")
print("\n== combined ==")
print(out.to_string(index=False))
