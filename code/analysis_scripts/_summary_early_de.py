import pandas as pd, numpy as np, os
root=r"D:\桌面\single cell\analysis_output"
d=pd.read_csv(os.path.join(root,"PH_early_treated_vs_blank_DE.csv"))
print(d.columns.tolist())
sig=d[(d.pvals_adj<0.05)&(d.logfoldchanges.abs()>0.5)].copy()
print("sig",len(sig),"up",(sig.logfoldchanges>0).sum(),"down",(sig.logfoldchanges<0).sum())
for title,dd in [("UP",sig.sort_values("scores",ascending=False).head(20)),("DOWN",sig.sort_values("scores").head(20))]:
 print("\n",title)
 print(dd[["gene","logfoldchanges","pvals_adj","pct_nz_group","scores"]].to_string(index=False))
