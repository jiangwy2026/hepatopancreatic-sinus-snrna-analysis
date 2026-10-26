import pandas as pd, os
f=r"D:\桌面\single cell\文章所需\Supplemental Table 5.xlsx"
d=pd.read_excel(f,sheet_name="shrimp-human orthologs")
print(d.shape); print(list(d.columns)); print(d.head(5).to_string())
for c in d.columns:
 print(c, d[c].astype(str).head(2).tolist())
