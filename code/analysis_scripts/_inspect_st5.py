import pandas as pd, os
f=r"D:\桌面\single cell\文章所需\Supplemental Table 5.xlsx"
for s in ["tpm>100","shrimp-human orthologs"]:
 d=pd.read_excel(f,sheet_name=s)
 print(s,d.shape); print(list(d.columns))
 print(d.head(3).to_string())
