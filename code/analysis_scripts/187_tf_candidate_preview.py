# -*- coding: utf-8 -*-
"""187_tf_candidate_preview.py - candidate TF evidence for PH proliferation."""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
OUT = r"D:\桌面\single cell\analysis_output"
PREV = os.path.join(OUT, "figures_pub", "sc_summary_new_names", "individual_previews")
os.makedirs(PREV, exist_ok=True)

E = pd.read_csv(os.path.join(OUT, "PH_TF_selected_evidence.csv"))
# keep the interpretable candidate module, ordered by target-program rank
keep = ["br", "Trl", "Cfp1", "Taf7", "Kdm4B", "Scm", "foxo", "Stat5b"]
E = E[E.TF.isin(keep)].copy()
E["order"] = E.TF.map({g:i for i,g in enumerate(keep)})
E = E.sort_values("order").reset_index(drop=True)

fig = plt.figure(figsize=(7.2, 3.25), facecolor="white")
# panel a: target-program rank enrichment
axa = fig.add_axes([0.07, 0.23, 0.36, 0.62])
yy = np.arange(len(E))
vals = E.mean_rank_pct.fillna(0.5).values
axa.barh(yy, vals, color="#4C78A8", height=0.58)
axa.axvline(0.5, color="#777777", lw=0.5, ls=(0, (2, 2)))
axa.set_yticks(yy); axa.set_yticklabels(E.TF, fontsize=5.6)
axa.invert_yaxis()
axa.set_xlim(0.35, 0.60); axa.set_xlabel("Mean target rank percentile", fontsize=5.7)
axa.set_title("TF target enrichment", fontsize=6.8, pad=4)
axa.tick_params(axis="x", labelsize=5.0, length=1.2); axa.tick_params(axis="y", length=0)
for i, r in E.iterrows():
    if np.isfinite(r.top10pct_frac):
        axa.text(r.mean_rank_pct + 0.006, i, "%.0f%%" % (100 * r.top10pct_frac), va="center", fontsize=4.9)
for sp in ("top", "right"): axa.spines[sp].set_visible(False)

# panel b: evidence heatmap
cols = ["early_logFC_proxy", "rho_G2M", "rho_ccnd2", "rho_pseudotime", "mean_rank_pct", "G2M_high_minus_low", "neglog10p"]
labels = ["Treated/control\nlog2FC", "rho G2M", "rho ccnd2", "rho\npseudotime", "Target rank\npercentile", "Target G2M\nhigh-low", "Target\n-enrichment"]
Z = E[cols].astype(float).copy()
Z["rho_pseudotime"] = -Z["rho_pseudotime"]  # higher = more progenitor-associated
# z-score each column for visual comparison
for c in cols:
    v = Z[c]
    if v.notna().sum() > 1:
        Z[c] = (v - v.mean(skipna=True)) / (v.std(skipna=True) + 1e-9)
axb = fig.add_axes([0.51, 0.23, 0.42, 0.62])
im = axb.imshow(Z.values, aspect="auto", cmap="RdBu_r", vmin=-1.8, vmax=1.8)
axb.set_xticks(range(len(cols))); axb.set_xticklabels(labels, fontsize=4.9, rotation=42, ha="right")
axb.set_yticks(range(len(E))); axb.set_yticklabels([""] * len(E))
axb.tick_params(length=0)
axb.set_title("Candidate TF evidence", fontsize=6.8, pad=4)
for i in range(len(E)):
    for j, c in enumerate(cols):
        v = E.loc[i, c] if c in E.columns else np.nan
        if c == "rho_pseudotime": v = -v
        txt = "NA" if pd.isna(v) else ("%.2f" % v if abs(v) < 10 else "%.1f" % v)
        axb.text(j, i, txt, ha="center", va="center", fontsize=4.5,
                 color="white" if abs(Z.iloc[i, j]) > 1.0 else "#222222")
cb = fig.colorbar(im, ax=axb, fraction=0.035, pad=0.02)
cb.set_label("Row-scaled evidence", fontsize=5.0); cb.ax.tick_params(labelsize=4.7, length=0.8); cb.outline.set_visible(False)
fig.text(0.5, 0.04,
         "No single dominant TF; br/Trl/Cfp1/Taf7/Kdm4B form a candidate chromatin/transcriptional activation module.",
         ha="center", va="center", fontsize=5.0, color="#333333")
p = os.path.join(PREV, "05_TF_candidates_preview.png")
fig.savefig(p, dpi=400, bbox_inches="tight", pad_inches=0.04, facecolor="white")
print("saved", p)
plt.close(fig)
