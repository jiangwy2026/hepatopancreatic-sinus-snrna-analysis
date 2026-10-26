import pandas as pd, os
files=[r"D:\桌面\single cell\文章所需\Table1.xlsx", r"D:\桌面\single cell\文章所需\图4数据\时序分析_119个蛋白集合汇总.xlsx", r"D:\桌面\single cell\文章所需\图3数据\正负离子差异代谢物汇总.xlsx", r"D:\桌面\single cell\文章所需\图3数据\代谢物汇总.pzfx"]
for f in files:
 print("\n==",f,os.path.exists(f))
 if f.endswith('.xlsx'):
  x=pd.ExcelFile(f); print("sheets",x.sheet_names)
  for s in x.sheet_names[:5]:
   d=pd.read_excel(f,sheet_name=s,nrows=5); print("--",s,d.shape); print(d.head().to_string())
