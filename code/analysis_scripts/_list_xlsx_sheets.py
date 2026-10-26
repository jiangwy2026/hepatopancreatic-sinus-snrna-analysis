import pandas as pd, os, glob
base=r"D:\桌面\single cell\文章所需"
files=glob.glob(os.path.join(base,"*.xlsx"))
for f in files:
 try:
  x=pd.ExcelFile(f)
  print(os.path.basename(f),x.sheet_names)
 except Exception as e: print("ERR",f,e)
