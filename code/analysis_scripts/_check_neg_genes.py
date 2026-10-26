import pandas as pd, numpy as np, os
D=pd.read_csv(r"D:\桌面\single cell\analysis_output\PH_early_treated_vs_blank_DE.csv")
D=D[np.isfinite(D.logfoldchanges)].copy()
sig=D[(D.pvals_adj<0.05)&(D.logfoldchanges<-0.5)]
for g in ["THBS4","PEN-3","LOC113801825","LOC113810905","LOC113830566","CYPA"]:
 print("\n",g)
 print(sig[sig.gene==g][["names","gene","logfoldchanges","pvals_adj","pct_nz_group","scores"]].round(3).to_string(index=False))
