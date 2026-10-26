import pandas as pd, os
O=r"D:\桌面\single cell\analysis_output"
S=pd.read_csv(os.path.join(O,"Fig6_marker_sensitivity.csv"))
for (p,o),d in S.dropna(subset=["best_group"]).groupby(["program","organ"]):
 g=d.sort_values("topN").best_group.tolist()
 if len(set(g))==1: print("STABLE",p,o,g[0])
 else: print("UNSTABLE",p,o,g)
