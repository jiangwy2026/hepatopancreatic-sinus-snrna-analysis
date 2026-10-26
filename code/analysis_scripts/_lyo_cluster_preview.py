import pandas as pd
import os
f=r"D:\桌面\single cell\DeGene.list.xls"
d=pd.read_excel(f,sheet_name="Top50.DeGene.list")
for c in sorted(d.Target_Cluster.unique()):
 x=d[d.Target_Cluster==c].sort_values("Log2FC",ascending=False).head(10)
 print("\nCluster",c,"n",len(d[d.Target_Cluster==c]))
 print(", ".join(x["Gene Name"].astype(str).tolist()))
