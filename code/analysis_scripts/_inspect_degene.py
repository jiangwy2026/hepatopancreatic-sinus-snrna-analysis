import pandas as pd, os
f=r"D:\桌面\single cell\DeGene.list.xls"
x=pd.ExcelFile(f)
print("sheets",x.sheet_names)
for s in x.sheet_names:
 d=pd.read_excel(f,sheet_name=s)
 print("\n==",s,d.shape); print(d.head(12).to_string())
 print("columns",d.columns.tolist())
