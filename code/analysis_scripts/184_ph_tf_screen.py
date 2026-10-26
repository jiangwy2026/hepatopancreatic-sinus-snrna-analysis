# -*- coding: utf-8 -*-
"""184_ph_tf_screen.py
Screen TFs for the PH proliferation response using:
- treated vs control in the early PH pseudotime bin;
- correlation with G2M score and ccnd2;
- negative correlation with differentiation pseudotime;
- TF-target co-expression programme score from SCENIC GRNBoost2 adjacencies.
This is candidate prioritisation, not causal proof.
"""
import os, sys, json
import numpy as np, pandas as pd, scanpy as sc
from scipy.stats import spearmanr, zscore

sc.settings.verbosity = 0
ROOT = r"D:\桌面\single cell"
OUT = os.path.join(ROOT, "analysis_output")

# --- annotations and PH object ---
ann = pd.read_csv(os.path.join(ROOT, "all.genes.marker.stat.xls.xls"), sep="\t", low_memory=False, encoding="gb18030")
ann.columns = [str(c).strip() for c in ann.columns]
tf = ann[ann["TF_family"].notna() & (ann["TF_family"].astype(str).str.strip() != "-") & (ann["TF_family"].astype(str).str.strip() != "")].copy()
tf["GeneID"] = tf["GeneID"].astype(str); tf["Symbol"] = tf["Symbol"].astype(str); tf["TF_family"] = tf["TF_family"].astype(str)
print("annotated TFs:", len(tf))

A = sc.read_h5ad(os.path.join(OUT, "adata_C2_within.h5ad"))
A.X = A.layers["counts"].copy(); sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
gid2i = {str(g): i for i, g in enumerate(A.var_names.astype(str))}
sym2gid = {}
for gid, s in zip(A.var_names.astype(str), A.var["sym"].astype(str)):
    sym2gid.setdefault(s, gid)
tf = tf[tf["GeneID"].isin(gid2i)].copy()
print("TFs in PH object:", len(tf))

# --- early DE table ---
DE = pd.read_csv(os.path.join(OUT, "PH_early_treated_vs_blank_DE.csv"))
DE["names"] = DE["names"].astype(str)
de_map = DE.set_index("names")[["logfoldchanges", "pvals_adj", "scores", "pct_nz_group"]]

# --- cell-level metrics ---
A.obs["pt_bin3"] = pd.qcut(A.obs["dpt_pseudotime"].astype(float), 3, labels=["P1", "P2", "P3"])
early = (A.obs["pt_bin3"] == "P1").values
group = A.obs["group"].astype(str).values
x = np.asarray(A.X.toarray(), dtype=np.float32)
cc_gid = sym2gid.get("ccnd2")
cc = x[:, gid2i[cc_gid]] if cc_gid in gid2i else np.full(A.n_obs, np.nan)
g2m = A.obs["G2M_score"].astype(float).values
pt = A.obs["dpt_pseudotime"].astype(float).values

# --- target programmes from GRNBoost2 adjacencies ---
adj = pd.read_csv(os.path.join(OUT, "scenic_adjacencies.csv"))
adj.columns = ["TF", "target", "importance"]
# map target symbols to gene IDs; duplicates keep first occurrence
adj["target_gid"] = adj["target"].astype(str).map(sym2gid)
adj = adj[adj["target_gid"].notna()].copy()
# z-score each gene across PH cells for programme scoring
Xz = (x - x.mean(0, keepdims=True)) / (x.std(0, keepdims=True) + 1e-6)
target_rows = []
for name, d in adj.groupby("TF"):
    d = d.sort_values("importance", ascending=False).drop_duplicates("target_gid").head(30)
    if len(d) < 10: continue
    idx = [gid2i[g] for g in d["target_gid"] if g in gid2i]
    if len(idx) < 10: continue
    score = Xz[:, idx].mean(1)
    r_g2m = spearmanr(score, g2m, nan_policy="omit").correlation
    r_cc = spearmanr(score, cc, nan_policy="omit").correlation
    r_pt = spearmanr(score, pt, nan_policy="omit").correlation
    target_rows.append({"TF": name, "n_targets": len(idx),
                        "target_rho_G2M": r_g2m, "target_rho_ccnd2": r_cc,
                        "target_rho_pseudotime": r_pt,
                        "target_mean_early_blank": float(score[early & (group == "blank")].mean()),
                        "target_mean_early_treated": float(score[early & (group == "treated")].mean())})
T = pd.DataFrame(target_rows)
T["target_diff_early"] = T["target_mean_early_treated"] - T["target_mean_early_blank"]
print("TF target programmes:", len(T))

# --- per-TF expression metrics ---
rows = []
for r in tf.itertuples():
    gid = str(r.GeneID); sym = str(r.Symbol); fam = str(r.TF_family)
    if gid not in gid2i: continue
    v = x[:, gid2i[gid]]
    if (v > 0).mean() < 0.05: continue
    rho_g2m = spearmanr(v, g2m, nan_policy="omit").correlation
    rho_cc = spearmanr(v, cc, nan_policy="omit").correlation
    rho_pt = spearmanr(v, pt, nan_policy="omit").correlation
    rec = {"TF": sym, "GeneID": gid, "TF_family": fam,
           "pct_detected": float((v > 0).mean() * 100),
           "mean_logexpr": float(v.mean()),
           "rho_G2M": rho_g2m, "rho_ccnd2": rho_cc, "rho_pseudotime": rho_pt,
           "early_blank_mean": float(v[early & (group == "blank")].mean()),
           "early_treated_mean": float(v[early & (group == "treated")].mean())}
    if gid in de_map.index:
        rec.update({"early_logFC": float(de_map.loc[gid, "logfoldchanges"]),
                    "early_pvals_adj": float(de_map.loc[gid, "pvals_adj"]),
                    "early_scores": float(de_map.loc[gid, "scores"]),
                    "early_pct_nz": float(de_map.loc[gid, "pct_nz_group"])})
    else:
        rec.update({"early_logFC": np.nan, "early_pvals_adj": np.nan, "early_scores": np.nan, "early_pct_nz": np.nan})
    rows.append(rec)
R = pd.DataFrame(rows)
# target programme merge by TF symbol
T = T.rename(columns={"TF": "TF"})
R = R.merge(T, on="TF", how="left")

# --- composite evidence score ---
def z(col):
    v = pd.to_numeric(R[col], errors="coerce")
    return (v - v.mean(skipna=True)) / v.std(skipna=True)
R["z_early_logFC"] = z("early_logFC")
R["z_rho_G2M"] = z("rho_G2M")
R["z_rho_ccnd2"] = z("rho_ccnd2")
R["z_neg_rho_pseudotime"] = z("rho_pseudotime") * -1
R["z_target_G2M"] = z("target_rho_G2M")
R["z_target_ccnd2"] = z("target_rho_ccnd2")
R["evidence_score"] = R[["z_early_logFC", "z_rho_G2M", "z_rho_ccnd2", "z_neg_rho_pseudotime", "z_target_G2M", "z_target_ccnd2"]].sum(axis=1, min_count=3)
R = R.sort_values("evidence_score", ascending=False).reset_index(drop=True)
R.to_csv(os.path.join(OUT, "PH_TF_candidates.csv"), index=False, encoding="utf-8-sig")
T.to_csv(os.path.join(OUT, "PH_TF_target_programs.csv"), index=False, encoding="utf-8-sig")
print("\nTop 30 TF candidates:")
cols = ["TF", "TF_family", "early_logFC", "early_pvals_adj", "rho_G2M", "rho_ccnd2", "rho_pseudotime", "target_rho_G2M", "target_rho_ccnd2", "evidence_score"]
print(R[cols].head(30).round(3).to_string(index=False))
json.dump({"n_annotated_TFs": int(len(tf)), "n_screened_TFs": int(len(R)),
           "n_target_programmes": int(len(T))}, open(os.path.join(OUT, "PH_TF_analysis_summary.json"), "w"), indent=1)
print("DONE")

