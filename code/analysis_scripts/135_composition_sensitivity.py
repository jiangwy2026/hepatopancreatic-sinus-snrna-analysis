import io, sys, os
import numpy as np, pandas as pd, anndata as ad
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
OUT = r"D:\桌面\single cell\analysis_output"
ORD = ["PH","GH","MH1","MH2","Epi","F","R","B","NEU"]
B = ad.read_h5ad(os.path.join(OUT, "adata_pertypeQC_harmony_tsne.h5ad"), backed="r")
obs = B.obs.copy()
# 方案1：厂商 QC（9,658）
V = obs[obs["orig_kept"].astype(bool)]
# 方案2：按类型阈值（18,638）
P = obs
def comp(d):
    t = pd.crosstab(d["ct"], d["sample"])
    return (t / t.sum(axis=0) * 100).reindex(ORD)
cv, cp = comp(V), comp(P)
tab = pd.DataFrame({"厂商QC_blank%": cv["blank"], "厂商QC_treated%": cv["treated"],
                    "按类型QC_blank%": cp["blank"], "按类型QC_treated%": cp["treated"]})
tab["倍数_厂商"] = (tab["厂商QC_treated%"] / tab["厂商QC_blank%"]).round(2)
tab["倍数_按类型"] = (tab["按类型QC_treated%"] / tab["按类型QC_blank%"]).round(2)
print("=== 两种 QC 下的亚群组成（% 为各样本内占比）===")
print(tab.round(2).to_string())
print()
print("样本细胞数：厂商 QC  blank %d / treated %d" % ((V["sample"]=="blank").sum(), (V["sample"]=="treated").sum()))
print("              按类型  blank %d / treated %d" % ((P["sample"]=="blank").sum(), (P["sample"]=="treated").sum()))
print()
print("=== 各样本细胞数按类型（按类型 QC）===")
print(pd.crosstab(P["ct"], P["sample"]).reindex(ORD).to_string())
print()
print("=== 各样本细胞数按类型（厂商 QC）===")
print(pd.crosstab(V["ct"], V["sample"]).reindex(ORD).to_string())
tab.to_csv(os.path.join(OUT, "composition_qc_sensitivity.csv"))
