import pandas as pd, os, glob
base=r"D:\桌面\single cell\文章所需"
for f in [os.path.join(base,"Supplemental Table 3.xlsx"),os.path.join(base,"Supplemental Table 5.xlsx"),os.path.join(base,"Supplemental Table 6.xlsx"),os.path.join(base,"Table1.xlsx")]:
 print("\n==",os.path.basename(f))
 x=pd.ExcelFile(f)
 for s in x.sheet_names:
  d=pd.read_excel(f,sheet_name=s,nrows=8)
  print("--",s,d.shape); print(d.head(8).to_string())
