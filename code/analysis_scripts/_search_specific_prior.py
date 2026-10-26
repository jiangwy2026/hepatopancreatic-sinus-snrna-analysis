import pandas as pd
a=pd.read_csv(r"D:\桌面\single cell\all.genes.marker.stat.xls.xls",sep="\t",encoding="gb18030",low_memory=False)
a.columns=[str(c).strip() for c in a.columns]
syms=["ARHGAP12","RPS3A","RPS9","EEF1A1","EEF1A","FLNA","FLN","PARVB","PARVA","SLC6A12","SLC6A12","BGT1","CHDH","BHMT","NFAT5","SOD2","SOD1","FTH1","GS2","GLUL","CHT","BBOX1","ALDH9A1","EIF4A2","EIF4E","EIF4H","EIF4G1","DYNLL1","DYNLL2","CTP","LOC113804001"]
for s in syms:
 m=a[a.Symbol.astype(str).str.lower()==s.lower()]
 if len(m):
  print("\n",s); print(m[["GeneID","Symbol","Description"]].to_string(index=False))
 else:
  m=a[a.Symbol.astype(str).str.lower().str.contains(s.lower(),regex=False)]
  if len(m): print("\n~",s); print(m[["GeneID","Symbol","Description"]].head(8).to_string(index=False))
