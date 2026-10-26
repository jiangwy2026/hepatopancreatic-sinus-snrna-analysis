import pandas as pd, os
t=pd.read_excel(r"D:\桌面\single cell\文章所需\Table1.xlsx",header=0)
t.to_csv(r"D:\桌面\single cell\analysis_output\Table1_rows.csv",index=False,encoding="utf-8-sig")
print(t.shape)
