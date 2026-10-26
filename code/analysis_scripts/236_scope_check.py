import pandas as pd, sys
sys.stdout.reconfigure(encoding="utf-8")
P1 = r"D:\桌面\single cell\文章所需\图5数据\肝胰腺比对材料.xlsx"
P2 = r"D:\桌面\single cell\文章所需\图5数据\淋巴比对材料.xlsx"
t = pd.read_excel(P1, sheet_name="tpm>100")
o = pd.read_excel(P1, sheet_name="shrimp-human orthologs")
l = pd.read_excel(P2, sheet_name="Sheet1")
lo = pd.read_excel(P2, sheet_name="shrimp-human orthologs")
print("tpm>100", t.shape, "ids", t.iloc[:,0].nunique())
print("lymphoid", l.shape)
tg = {
 "CXCL12":None,"KITLG":None,"LYVE1":None,"PTN":None,"DAB2":"113808873","KIT":"113804259",
 "CXCR4":"113819525","JAG1":"113821858","DLL1":"113815864","NOTCH1":"113802681",
 "NOTCH2":"113802400","WNT5A":"113806005","FZD5":"113826277","FZD4":"113826275",
 "ANGPT2":"113804527","TEK/TIE1":"113823485","THBS4":"113823118","LAMA1":"113802906",
 "ITGAV":"113821319","ITGB1":"113830322","ITGA4":"113827435","MRC1":"113822862",
 "CSF1R":"113823485","POSTN":"113805211","PROM1":"113802277","VWF":"113824171",
 "CDH5":"113802417","VEGFA":"113804841"}
ids = set(str(x).replace("ncbi_","") for x in t.iloc[:,0].astype(str))
lids = set(str(x).replace("ncbi_","") for x in l.iloc[:,0].astype(str))
for g, gid in tg.items():
    print(f"{g:9s} gid={str(gid):10s} in_tpm100={gid in ids if gid else 'NA'} in_lymphoid={gid in lids if gid else 'NA'}")
