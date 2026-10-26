import pandas as pd, numpy as np, os, re
O=r"D:\桌面\single cell\analysis_output"
D=pd.read_csv(os.path.join(O,"PH_early_treated_vs_blank_DE.csv"))
D=D[np.isfinite(D.logfoldchanges)].copy()
pat=re.compile(r'(thbs|pen-3|crustin|wap|pten|p53|tp53|p21|p27|cdkn|foxo|smad|tgfb|bmp|dpp|noggin|gremlin|chordin|dkk|sfrp|wifs|axin|apc|gsk|hippo|hpo|wts|lats|yorkie|yki|merlin|nf2|kibra|rassf|tob|ski|sno|rb1|retinoblastoma)',re.I)
M=D[D.gene.astype(str).str.contains(pat,na=False)]
print(M[["gene","logfoldchanges","pvals_adj","scores"]].sort_values("logfoldchanges").to_string(index=False))
