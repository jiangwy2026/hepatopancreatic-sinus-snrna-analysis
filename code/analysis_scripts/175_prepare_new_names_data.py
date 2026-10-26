# -*- coding: utf-8 -*-
"""175_prepare_new_names_data.py
Recompute the data needed for the new-name composition, velocity summary and
PH early-response mechanism panels.
"""
import os, sys, json
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.stats import chi2_contingency, mannwhitneyu

sc.settings.verbosity = 0
OUT = r"D:\桌面\single cell\analysis_output"
GROUPS = {"0": "TS", "1": "GH", "2": "PH", "3": "MH2", "4": "PM-like",
          "5": "F", "6": "R", "7": "B", "8": "MyoEpi", "9": "M"}
OLD2NEW = {"PH": "PH", "GH": "GH", "MH1": "TS", "MH2": "MH2"}

# ---------------- composition ----------------
A = ad.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"))
ct = A.obs["vendor_cluster"].astype(str).map(GROUPS)
sample = A.obs["sample"].astype(str).values
rows = []
for g in sorted(GROUPS.values(), key=lambda x: int((ct == x).sum()), reverse=True):
    nb = int(((sample == "blank") & (ct == g)).sum())
    nt = int(((sample == "treated") & (ct == g)).sum())
    pb = 100.0 * nb / int((sample == "blank").sum())
    pt = 100.0 * nt / int((sample == "treated").sum())
    rows.append({"group": g, "n_blank": nb, "n_treated": nt, "n_total": nb + nt,
                 "pct_blank": round(pb, 3), "pct_treated": round(pt, 3),
                 "fold_treated_vs_blank": round(pt / pb, 3) if pb > 0 else np.nan})
comp = pd.DataFrame(rows)
comp.to_csv(os.path.join(OUT, "composition_new_names_10groups.csv"), index=False, encoding="utf-8-sig")
tab = np.array([[int(((sample == "blank") & (ct == g)).sum()) for g in comp.group],
                [int(((sample == "treated") & (ct == g)).sum()) for g in comp.group]])
chi2, p, dof, _ = chi2_contingency(tab)
json.dump({"chi2": float(chi2), "p": float(p), "dof": int(dof),
           "n_blank": int((sample == "blank").sum()), "n_treated": int((sample == "treated").sum())},
          open(os.path.join(OUT, "composition_new_names_chi2.json"), "w"), indent=1)
print("composition:\n", comp.to_string(index=False))
print("chi2", chi2, "p", p, "dof", dof)

# ---------------- velocity ----------------
V = pd.read_csv(os.path.join(OUT, "velocity_group_flow_normalized.csv"), index_col=0)
V = V.rename(index=OLD2NEW, columns=OLD2NEW)
vel_order = ["PH", "GH", "TS", "MH2"]
V4 = V.loc[vel_order, vel_order].copy()
V4.to_csv(os.path.join(OUT, "velocity_flow_4groups_new_names.csv"), encoding="utf-8-sig")
Voff = V4.astype(float).copy()
_arr = Voff.to_numpy(copy=True)
np.fill_diagonal(_arr, np.nan)
Voff = pd.DataFrame(_arr, index=V4.index, columns=V4.columns)
Voff.to_csv(os.path.join(OUT, "velocity_flow_4groups_offdiag.csv"), encoding="utf-8-sig")
print("velocity flow:\n", V4.round(4).to_string())

VP = pd.read_csv(os.path.join(OUT, "velocity_per_group.csv"))
VP["group_new"] = VP["group"].map(OLD2NEW).fillna(VP["group"])
VP = VP[VP["group_new"].isin(vel_order)].copy()
VP.to_csv(os.path.join(OUT, "velocity_per_group_new_names.csv"), index=False, encoding="utf-8-sig")
print("velocity pseudotime:\n", VP[["group_new", "velocity_pseudotime", "velocity_confidence"]].to_string(index=False))

# ---------------- PH early response ----------------
A2 = ad.read_h5ad(os.path.join(OUT, "adata_C2_within.h5ad"))
A2.X = A2.layers["counts"].copy()
sc.pp.normalize_total(A2, target_sum=1e4)
sc.pp.log1p(A2)
A2.obs["pt_bin3"] = pd.qcut(A2.obs["dpt_pseudotime"].astype(float), 3, labels=["P1", "P2", "P3"])
early = A2[A2.obs["pt_bin3"] == "P1"].copy()
print("early PH:", early.n_obs, early.obs["group"].value_counts().to_dict())
sc.tl.rank_genes_groups(early, "group", groups=["treated"], reference="blank", method="wilcoxon", pts=True)
d = sc.get.rank_genes_groups_df(early, group="treated")
sym = dict(zip(early.var_names.astype(str), early.var["sym"].astype(str)))
d["sym"] = d["names"].astype(str).map(sym)
d["gene"] = d["sym"].where(d["sym"].notna() & (d["sym"] != "nan"), d["names"].astype(str))
d.to_csv(os.path.join(OUT, "PH_early_treated_vs_blank_DE.csv"), index=False, encoding="utf-8-sig")
print("early DE rows:", len(d))

sym2gid = {}
for gid, s in zip(A2.var_names.astype(str), A2.var["sym"].astype(str)):
    sym2gid.setdefault(s, gid)
stats = []
for symname, label in [("ccnd2", "ccnd2"), ("LOC113823945", "hemocyte transglutaminase (LOC113823945)"),
                       ("TGM1", "TGM1")]:
    gid = sym2gid.get(symname)
    if gid is None:
        continue
    v = early[:, gid].X.toarray().ravel().astype(float)
    for grp in ("blank", "treated"):
        m = (early.obs["group"].astype(str) == grp).values
        vv = v[m]
        stats.append({"gene": symname, "label": label, "group": grp, "n": int(m.sum()),
                      "mean_logexpr": float(vv.mean()), "median_logexpr": float(np.median(vv)),
                      "pct_detected": float((vv > 0).mean() * 100)})
    mb = (early.obs["group"].astype(str) == "blank").values
    mt = (early.obs["group"].astype(str) == "treated").values
    u, pv = mannwhitneyu(v[mt], v[mb], alternative="two-sided")
    stats.append({"gene": symname, "label": label, "group": "p_value", "n": int(early.n_obs),
                  "mean_logexpr": np.nan, "median_logexpr": np.nan, "pct_detected": np.nan,
                  "mannwhitney_p": float(pv)})
S = pd.DataFrame(stats)
S.to_csv(os.path.join(OUT, "PH_early_ccnd2_TGM_stats.csv"), index=False, encoding="utf-8-sig")
print("mechanism stats:\n", S.to_string(index=False))
print("DONE")



