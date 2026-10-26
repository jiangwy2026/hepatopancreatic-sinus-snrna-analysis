import pandas as pd, glob, os
for f in sorted(glob.glob(r"D:\桌面\single cell\analysis_output\elife_ref\*.csv")):
 print("\n==",os.path.basename(f))
 d=pd.read_csv(f)
 print(d.shape); print(d.columns.tolist()); print(d.head(5).to_string())
