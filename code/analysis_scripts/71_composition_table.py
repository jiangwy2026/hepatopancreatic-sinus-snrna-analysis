import sys, io, os
import numpy as np, pandas as pd, scanpy as sc
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 250); sc.settings.verbosity = 0
out = r"D:\桌面\single cell\analysis_output"
A = sc.read_h5ad(out + r"\adata_soupx.h5ad")
base = {"2":"PH","1":"GH","0":"MH1","3":"MH2","4":"Epi","5":"F","6":"R","7":"B","8":"NE","9":"Neu"}
A.obs["ct"] = A.obs.vendor_cluster.astype(str).map(base).replace({"NE":"NEU","Neu":"NEU"})
ORD = ["PH","GH","MH1","MH2","Epi","F","R","B","NEU"]
A.obs["ct"] = pd.Categorical(A.obs["ct"], categories=ORD)
t = pd.crosstab(A.obs["ct"], A.obs["group"])
t = t[["blank","treated"]]
pct = t / t.sum(0) * 100
pct["倍数"] = (t["treated"]/t["treated"].sum()) / (t["blank"]/t["blank"].sum())
print("=== 细胞数 ==="); print(t.to_string())
print("\n=== 百分比 (%) ==="); print(pct.round(2).to_string())
pct.to_csv(out + r"\composition_pct_by_group.csv")
# 卡方检验（整体组成是否不同）
from scipy.stats import chi2_contingency
chi2, p, dof, _ = chi2_contingency(t.values)
print("\n整体组成卡方检验: chi2=%.1f, dof=%d, p=%.3e" % (chi2, dof, p))
