import pandas as pd, glob, os
for f in sorted(glob.glob(r"D:\桌面\single cell\analysis_output\elife_ref\*.csv")):
 d=pd.read_csv(f)
 print(os.path.basename(f), sorted(d.Target_Cluster.unique().tolist()))
