import pandas as pd, os
for f in [r"D:\桌面\single cell\文章所需\图5数据\淋巴比对材料.xlsx",r"D:\桌面\single cell\文章所需\图5数据\肝血窦-淋巴基因汇总.xlsx",r"D:\桌面\single cell\文章所需\Supplemental Table 6.xlsx"]:
 print("\n==",f)
 x=pd.ExcelFile(f); print(x.sheet_names)
 for s in x.sheet_names[:6]:
  d=pd.read_excel(f,sheet_name=s,nrows=8); print("--",s,d.shape); print(d.head(8).to_string()); print(d.columns.tolist())
