# -*- coding: utf-8 -*-
"""Cross-organ cell-cycle-adjusted progenitor potential.

Primary contract
-----------------
Core conclusion: determine whether the hepatic-sinus PH population retains a
higher non-cycling progenitor-identity program than proliferative prohemocytes
from circulating hemolymph and the PH-like R6 population of the hepatopancreas.

Evidence chain:
1. UCell scores for the eLife PH1 cell-cycle program and PH2 identity program.
2. Within-dataset percentile normalization to remove platform/scale differences.
3. Background-only, cell-cycle-binned prediction of the PH2 program.
4. Standardized residual ("progenitor potential beyond cell cycle").
5. Cycle-matched permutation and stratified bootstrap intervals.

The analysis is cell-level and exploratory; cell counts are not biological
replicates. LyO has no cell-level matrix in the workspace and is therefore not
included as a fourth organ in the adjusted comparison.
"""

from pathlib import Path
import json
import numpy as np
import pandas as pd
import anndata as ad
from scipy.stats import rankdata, pearsonr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import Patch

ROOT = Path(r"D:\桌面\single cell")
OUT = ROOT / "analysis_output" / "cross_organ_progenitor"
OUT.mkdir(parents=True, exist_ok=True)
HS_PATH = ROOT / "analysis_output" / "step3" / "step3_lyo_scores.h5ad"
ELIFE_PATH = ROOT / "analysis_output" / "elife_adata.h5ad"
GEO_RAW_PATH = ROOT / "analysis_output" / "geo_GSE306674_raw.h5ad"
GEO_PROC_PATH = ROOT / "analysis_output" / "geo_GSE306674_processed.h5ad"
PH1_PATH = ROOT / "analysis_output" / "elife_ref" / "elife-80127-supp3-v2.csv"
PH2_PATH = ROOT / "analysis_output" / "elife_ref" / "elife-80127-supp6-v2.csv"

RNG = np.random.default_rng(20261021)
assert fm.findfont("Times New Roman", fallback_to_default=False)
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman"],
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "axes.linewidth": 0.65,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "savefig.facecolor": "white",
})


def norm_id(x):
    return str(x).replace("ncbi_", "")


def load_gene_ids(path, top_n=50):
    d = pd.read_csv(path)
    d = d.sort_values("Log2FC", ascending=False).head(top_n)
    return list(dict.fromkeys(norm_id(x) for x in d["Gene ID"].astype(str)))


def map_ids(gene_ids, var_names):
    lookup = {}
    for i, x in enumerate(map(norm_id, var_names)):
        lookup.setdefault(x, []).append(i)
    mapped, labels, missing = [], [], []
    for x in gene_ids:
        hit = lookup.get(x, [])
        if len(hit) == 1:
            mapped.append(hit[0]); labels.append(x)
        else:
            missing.append(x)
    return mapped, labels, missing


def ucell_scores(layer, row_positions, gene_idx, n_vars, chunk=500):
    row_positions = np.asarray(row_positions)
    out = np.zeros(len(row_positions), dtype=np.float32)
    k = len(gene_idx)
    if k == 0:
        raise ValueError("Empty gene set")
    denom = k * (n_vars - k)
    for st in range(0, len(row_positions), chunk):
        en = min(len(row_positions), st + chunk)
        pos = row_positions[st:en]
        X = layer[pos, :]
        if hasattr(X, "toarray"):
            X = np.asarray(X.toarray(), dtype=np.float32)
        else:
            X = np.asarray(X, dtype=np.float32)
        order = np.argsort(-X, axis=1, kind="stable")
        ranks = np.empty_like(order)
        ranks[np.arange(en - st)[:, None], order] = np.arange(1, n_vars + 1, dtype=np.int32)
        S = ranks[:, gene_idx].sum(axis=1)
        out[st:en] = 1.0 - (S - k * (k + 1) / 2.0) / denom
    return out


def pct_rank(x):
    return (rankdata(x, method="average") - 0.5) / len(x)


def fit_cycle_trend(cycle_pct, prog_pct, background, n_bins=20):
    edges = np.linspace(0, 1, n_bins + 1)
    centers = (edges[:-1] + edges[1:]) / 2
    ids = np.clip(np.digitize(cycle_pct, edges) - 1, 0, n_bins - 1)
    y = np.full(n_bins, np.nan, dtype=float)
    for b in range(n_bins):
        m = background & (ids == b)
        if m.sum() >= 40:
            y[b] = np.median(prog_pct[m])
    good = np.isfinite(y)
    if good.sum() < 3:
        raise ValueError("Too few populated cycle bins")
    y = np.interp(np.arange(n_bins), np.flatnonzero(good), y[good])
    pred = np.interp(cycle_pct, centers, y)
    return ids, centers, y, pred


def stratified_permutation(resid, target, cycle_pct, n_perm=10000, n_bins=10):
    edges = np.linspace(0, 1, n_bins + 1)
    ids = np.clip(np.digitize(cycle_pct, edges) - 1, 0, n_bins - 1)
    observed = float(resid[target].mean())
    null = np.empty(n_perm, dtype=np.float64)
    rng = np.random.default_rng(20261021)
    for p in range(n_perm):
        pseudo = np.zeros(len(resid), dtype=bool)
        for b in range(n_bins):
            idx = np.flatnonzero(ids == b)
            k = int(target[idx].sum())
            if k:
                pseudo[rng.choice(idx, size=k, replace=False)] = True
        null[p] = resid[pseudo].mean()
    pvalue = (1.0 + np.sum(null >= observed)) / (1.0 + n_perm)
    return observed, pvalue, null


def stratified_bootstrap(resid, target, cycle_pct, n_boot=5000, n_bins=10):
    edges = np.linspace(0, 1, n_bins + 1)
    ids = np.clip(np.digitize(cycle_pct, edges) - 1, 0, n_bins - 1)
    rng = np.random.default_rng(20261022)
    vals = np.empty(n_boot, dtype=np.float64)
    groups = []
    for b in range(n_bins):
        v = resid[target & (ids == b)]
        if len(v):
            groups.append(v)
    for i in range(n_boot):
        draws = [rng.choice(v, size=len(v), replace=True) for v in groups]
        vals[i] = np.concatenate(draws).mean()
    return np.quantile(vals, [0.025, 0.975])


def analyze_hs():
    a = ad.read_h5ad(HS_PATH, backed="r")
    cyc_idx, cyc_lab, cyc_missing = map_ids(CYCLE_IDS, a.var_names)
    prog_idx, prog_lab, prog_missing = map_ids(PROG_IDS, a.var_names)
    pos = np.arange(a.n_obs)
    cyc = ucell_scores(a.layers["counts"], pos, cyc_idx, a.n_vars)
    prog = ucell_scores(a.layers["counts"], pos, prog_idx, a.n_vars)
    obs = a.obs.copy()
    a.file.close()
    target = obs["cell_type"].astype(str).eq("PH").to_numpy()
    eligible = np.ones(a.n_obs, dtype=bool)
    return {
        "name": "Hs PH",
        "dataset": "Hepatic sinus",
        "group": "PH",
        "cycle": cyc,
        "prog": prog,
        "target": target,
        "eligible": eligible,
        "obs": obs,
        "cycle_mapped": len(cyc_idx),
        "prog_mapped": len(prog_idx),
        "cycle_missing": cyc_missing,
        "prog_missing": prog_missing,
    }


def analyze_elife():
    a = ad.read_h5ad(ELIFE_PATH, backed="r")
    cyc_idx, cyc_lab, cyc_missing = map_ids(CYCLE_IDS, a.var_names)
    prog_idx, prog_lab, prog_missing = map_ids(PROG_IDS, a.var_names)
    pos = np.arange(a.n_obs)
    cyc = ucell_scores(a.layers["counts"], pos, cyc_idx, a.n_vars)
    prog = ucell_scores(a.layers["counts"], pos, prog_idx, a.n_vars)
    obs = a.obs.copy()
    a.file.close()
    cl = obs["elife_cluster"].astype(object).where(obs["elife_cluster"].notna(), "").astype(str)
    eligible = (cl.ne("") & cl.ne("nan")).to_numpy()
    target = cl.isin(["8", "11"]).to_numpy()
    return {
        "name": "Blood PH",
        "dataset": "Circulating hemolymph",
        "group": "PH1-like + PH2-like (clusters 8,11)",
        "cycle": cyc,
        "prog": prog,
        "target": target,
        "eligible": eligible,
        "obs": obs,
        "cycle_mapped": len(cyc_idx),
        "prog_mapped": len(prog_idx),
        "cycle_missing": cyc_missing,
        "prog_missing": prog_missing,
    }


def analyze_geo():
    raw = ad.read_h5ad(GEO_RAW_PATH, backed="r")
    pro = ad.read_h5ad(GEO_PROC_PATH, backed="r")
    cyc_idx, cyc_lab, cyc_missing = map_ids(CYCLE_IDS, raw.var_names)
    prog_idx, prog_lab, prog_missing = map_ids(PROG_IDS, raw.var_names)
    pos = np.array([raw.obs_names.get_loc(b) for b in pro.obs_names], dtype=np.int64)
    cyc = ucell_scores(raw.layers["counts"], pos, cyc_idx, raw.n_vars)
    prog = ucell_scores(raw.layers["counts"], pos, prog_idx, raw.n_vars)
    obs = pro.obs.copy()
    raw.file.close(); pro.file.close()
    cl = obs["cl"].astype(str)
    target = cl.eq("6").to_numpy()
    eligible = np.ones(obs.shape[0], dtype=bool)
    return {
        "name": "Hepato R6",
        "dataset": "Hepatopancreas",
        "group": "PH-like R6 (cluster 6)",
        "cycle": cyc,
        "prog": prog,
        "target": target,
        "eligible": eligible,
        "obs": obs,
        "cycle_mapped": len(cyc_idx),
        "prog_mapped": len(prog_idx),
        "cycle_missing": cyc_missing,
        "prog_missing": prog_missing,
    }


CYCLE_IDS = load_gene_ids(PH1_PATH)
PROG_IDS = load_gene_ids(PH2_PATH)
print("Cycle genes", len(CYCLE_IDS), "Program genes", len(PROG_IDS))

results = []
percell = []
stratified_rows = []
for fn in (analyze_hs, analyze_elife, analyze_geo):
    d = fn()
    elig = d["eligible"]
    cyc = d["cycle"]
    prog = d["prog"]
    cyc_pct = pct_rank(cyc)
    prog_pct = pct_rank(prog)
    target = d["target"] & elig
    background = (~d["target"]) & elig
    bins, centers, trend, pred = fit_cycle_trend(cyc_pct, prog_pct, background)
    resid = prog_pct - pred
    resid_z = (resid - resid.mean()) / (resid.std(ddof=1) + 1e-12)
    obs_mean, perm_p, null = stratified_permutation(resid_z, target, cyc_pct)
    ci_low, ci_high = stratified_bootstrap(resid_z, target, cyc_pct)
    r_bg = pearsonr(cyc_pct[background], prog_pct[background])
    results.append({
        "target": d["name"], "dataset": d["dataset"], "target_definition": d["group"],
        "n_total": int(elig.sum()), "n_target": int(target.sum()),
        "cycle_genes_mapped": d["cycle_mapped"], "progenitor_genes_mapped": d["prog_mapped"],
        "mean_cycle_pct": float(cyc_pct[target].mean()),
        "mean_prog_pct_raw": float(prog_pct[target].mean()),
        "mean_adjusted_prog_z": float(obs_mean),
        "adjusted_ci_low": float(ci_low), "adjusted_ci_high": float(ci_high),
        "frac_adjusted_gt0": float((resid_z[target] > 0).mean()),
        "frac_adjusted_gt1": float((resid_z[target] > 1).mean()),
        "cycle_prog_r_background": float(r_bg.statistic),
        "cycle_prog_p_background": float(r_bg.pvalue),
        "permutation_p": float(perm_p),
        "permutation_n": 10000,
    })
    # per-cell values
    pp = pd.DataFrame({
        "cell": d["obs"].index.astype(str),
        "dataset": d["dataset"],
        "target": d["name"],
        "group": d["group"],
        "is_target": target,
        "cycle_ucell": cyc,
        "progenitor_ucell": prog,
        "cycle_pct": cyc_pct,
        "progenitor_pct": prog_pct,
        "predicted_progenitor_pct": pred,
        "adjusted_progenitor_z": resid_z,
    })
    percell.append(pp)
    # stratified target/background means
    edges = np.linspace(0, 1, 11)
    for b in range(10):
        inbin = (cyc_pct >= edges[b]) & (cyc_pct < edges[b + 1] if b < 9 else cyc_pct <= 1.0)
        for is_target, label in ((True, "target"), (False, "background")):
            m = inbin & (target if is_target else background)
            stratified_rows.append({
                "target": d["name"], "bin": b + 1,
                "bin_left": edges[b], "bin_right": edges[b + 1],
                "group": label, "n": int(m.sum()),
                "mean_cycle_pct": float(cyc_pct[m].mean()) if m.any() else np.nan,
                "mean_prog_pct": float(prog_pct[m].mean()) if m.any() else np.nan,
                "mean_adjusted_z": float(resid_z[m].mean()) if m.any() else np.nan,
            })
    print(d["name"], "n", target.sum(), "mapped", d["cycle_mapped"], d["prog_mapped"],
          "adjusted", round(obs_mean, 3), "CI", round(ci_low, 3), round(ci_high, 3), "p", perm_p)

summary = pd.DataFrame(results)
percell_df = pd.concat(percell, ignore_index=True)
strat = pd.DataFrame(stratified_rows)
summary.to_csv(OUT / "cross_organ_adjusted_summary.csv", index=False, encoding="utf-8-sig")
percell_df.to_csv(OUT / "cross_organ_percell_scores.csv.gz", index=False, compression="gzip")
strat.to_csv(OUT / "cross_organ_cycle_stratified_means.csv", index=False, encoding="utf-8-sig")
metadata = {
    "cycle_program": "eLife PH1 top50",
    "progenitor_program": "eLife PH2 top50 (non-cycling identity program)",
    "hemolymph_target": "eLife clusters 8 and 11 (PH1-like and PH2-like)",
    "hepatopancreas_target": "cluster 6 (R6, strongest PH label-transfer match)",
    "normalization": "within-dataset percentile",
    "adjustment": "20-bin background cycle trend; residual standardized within dataset",
    "test": "cycle-decile-stratified label permutation (10,000 permutations)",
    "interval": "cycle-decile-stratified bootstrap (5,000 resamples)",
}
(OUT / "analysis_metadata.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")

# Plot
labels = summary["target"].tolist()
colors = {"Hs PH": "#B2182B", "Blood PH": "#2166AC", "Hepato R6": "#4D9221"}
fig = plt.figure(figsize=(7.2, 5.1))
gs = fig.add_gridspec(2, 3, height_ratios=[1.15, 0.85], hspace=0.42, wspace=0.35)
# A: raw relation, three datasets
for j, name in enumerate(labels):
    ax = fig.add_subplot(gs[0, j])
    d = percell_df[percell_df.target == name]
    bg = d[~d.is_target]
    tg = d[d.is_target]
    ax.hexbin(bg.cycle_pct, bg.progenitor_pct, gridsize=35, cmap="Greys", mincnt=1, linewidths=0)
    ax.scatter(tg.cycle_pct, tg.progenitor_pct, s=4, alpha=0.25, color=colors[name], rasterized=True)
    row = summary[summary.target == name].iloc[0]
    # fitted trend from background
    xgrid = np.linspace(0, 1, 100)
    # reconstruct trend by bins as stored in stratified means
    sm = strat[(strat.target == name) & (strat.group == "background")]
    ax.plot(sm.mean_cycle_pct, sm.mean_prog_pct, color=colors[name], lw=1.2)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_title(f"{name}\n(n={int(row.n_target)})", pad=4)
    ax.set_xlabel("Cell-cycle percentile")
    if j == 0: ax.set_ylabel("Progenitor program percentile")
    ax.text(0.04, 0.95, "ABC"[j], transform=ax.transAxes, va="top", fontsize=10, fontweight="bold")
# B: forest plot adjusted potential
ax = fig.add_subplot(gs[1, :2])
y = np.arange(len(summary))[::-1]
for yi, (_, r) in zip(y, summary.iterrows()):
    ax.errorbar(r.mean_adjusted_prog_z, yi,
                xerr=[[r.mean_adjusted_prog_z - r.adjusted_ci_low], [r.adjusted_ci_high - r.mean_adjusted_prog_z]],
                fmt='o', ms=6, capsize=3, color=colors[r.target], ecolor="#555555", lw=1)
    ax.text(r.adjusted_ci_high + 0.06, yi, f"z={r.mean_adjusted_prog_z:.2f}, P={r.permutation_p:.3g}",
            va="center", fontsize=6.7)
ax.axvline(0, color="#777", lw=0.7, ls="--")
ax.set_yticks(y); ax.set_yticklabels(summary.target)
ax.set_xlabel("Adjusted progenitor potential (z)")
ax.set_title("Progenitor identity after cell-cycle adjustment", pad=5)
ax.text(-0.12, 1.02, "D", transform=ax.transAxes, fontsize=10, fontweight="bold")
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.set_xlim(min(-0.5, summary.adjusted_ci_low.min() - 0.2), summary.adjusted_ci_high.max() + 0.9)
# E: target versus background residual distributions
ax = fig.add_subplot(gs[1, 2])
positions = np.arange(len(labels))
for i, name in enumerate(labels):
    d = percell_df[percell_df.target == name]
    bg_vals = d.loc[~d.is_target, "adjusted_progenitor_z"].dropna().values
    tg_vals = d.loc[d.is_target, "adjusted_progenitor_z"].dropna().values
    vp1 = ax.violinplot([bg_vals], positions=[i - 0.18], widths=0.28, showextrema=False)
    vp2 = ax.violinplot([tg_vals], positions=[i + 0.18], widths=0.28, showextrema=False)
    for body in vp1["bodies"]:
        body.set_facecolor("#BDBDBD"); body.set_edgecolor("none"); body.set_alpha(0.7)
    for body in vp2["bodies"]:
        body.set_facecolor(colors[name]); body.set_edgecolor("none"); body.set_alpha(0.75)
    ax.plot([i - 0.18, i + 0.18], [np.median(bg_vals), np.median(tg_vals)],
            color="#333333", lw=0.7, zorder=3)
    ax.scatter([i - 0.18, i + 0.18], [np.median(bg_vals), np.median(tg_vals)],
               s=6, color=["#777777", colors[name]], zorder=4)
ax.axhline(0, color="#777777", lw=0.7, ls="--")
ax.set_xticks(positions); ax.set_xticklabels(labels, fontsize=6)
ax.set_ylabel("Adjusted progenitor potential (z)")
ax.set_title("Target versus background", pad=5)
ax.text(-0.18, 1.16, "E", transform=ax.transAxes, fontsize=10, fontweight="bold")
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.legend(handles=[Patch(facecolor="#BDBDBD", label="background"), Patch(facecolor="#B2182B", label="target (colored)")],
          frameon=False, fontsize=5.4, loc="upper left")
fig.savefig(OUT / "Fig_cross_organ_adjusted_progenitor.png", dpi=300, bbox_inches="tight")
fig.savefig(OUT / "Fig_cross_organ_adjusted_progenitor.svg", format="svg", bbox_inches="tight")
fig.savefig(OUT / "Fig_cross_organ_adjusted_progenitor.pdf", format="pdf", bbox_inches="tight")
fig.savefig(OUT / "Fig_cross_organ_adjusted_progenitor.tiff", format="tiff", dpi=1200, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
plt.close(fig)

print(summary.to_string(index=False))
print("saved", OUT)









