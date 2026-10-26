import pandas as pd, os
p=r"E:\workspace\scenic_work\expr.csv.gz"
print("exists",os.path.exists(p))
if os.path.exists(p):
 e=pd.read_csv(p,index_col=0,nrows=3)
 print(e.shape); print(e.columns[:10].tolist())
try:
 import pyarrow.feather as feather
 p2=r"D:\桌面\single cell\analysis_output\scenic_db\dm6_v10_clust.genes_vs_motifs.rankings.feather"
 t=feather.read_table(p2).to_pandas()
 print("ranking",t.shape); print(t["genes"].head(10).tolist() if "genes" in t.columns else t.head())
except Exception as ex:
 print("feather read failed",repr(ex))
