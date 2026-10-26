# -*- coding: utf-8 -*-
"""Preview a replacement for Fig. 2D: PH identity versus cell-cycle state."""

from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from scipy.stats import rankdata
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

ROOT = Path(r"D:\桌面\single cell")
OUT = ROOT / "analysis_output" / "figures_pub" / "2026-10-23_Fig2D_replacement"
OUT.mkdir(parents=True, exist_ok=True)
H5 = ROOT / "analysis_output" / "step3" / "step3_lyo_scores.h5ad"
PH1 = ROOT / "analysis_output" / "elife_ref" / "elife-80127-supp3-v2.csv"
PH2 = ROOT / "analysis_output" / "elife_ref" / "elife-80127-supp6-v2.csv"

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


def ids(path):
    d = pd.read_csv(path)
    return list(dict.fromkeys(d.sort_values("Log2FC", ascending=False).head(50)["Gene ID"].astype(str)))


def norm(x):
    return str(x).replace("ncbi_", "")


cycle_ids = [norm(x) for x in ids(PH1)]
prog_ids = [norm(x) for x in ids(PH2) if norm(x) not in set(cycle_ids)]

a = ad.read_h5ad(H5, backed="r")
var_ids = [norm(x) for x in a.var_names]
lookup = {v: i for i, v in enumerate(var_ids)}

cycle_idx = [lookup[x] for x in cycle_ids if x in lookup]
prog_idx = [lookup[x] for x in prog_ids if x in lookup]
print("mapped cycle", len(cycle_idx), "mapped identity", len(prog_idx))


def ucell(idx, chunk=500):
    n, N = a.n_obs, a.n_vars
    out = np.zeros(n, dtype=np.float32)
    k = len(idx)
    for st in range(0, n, chunk):
        en = min(n, st + chunk)
        X = a.layers["counts"][st:en, :]
        X = np.asarray(X.toarray(), dtype=np.float32) if hasattr(X, "toarray") else np.asarray(X, dtype=np.float32)
        order = np.argsort(-X, axis=1, kind="stable")
        ranks = np.empty_like(order)
        ranks[np.arange(en - st)[:, None], order] = np.arange(1, N + 1, dtype=np.int32)
        S = ranks[:, idx].sum(axis=1)
        out[st:en] = 1.0 - (S - k * (k + 1) / 2.0) / (k * (N - k))
    return out


cycle = ucell(cycle_idx)
prog = ucell(prog_idx)
obs = a.obs.copy()
a.file.close()

cycle_pct = (rankdata(cycle, method="average") - 0.5) / len(cycle)
prog_pct = (rankdata(prog, method="average") - 0.5) / len(prog)
ph = obs["cell_type"].astype(str).eq("PH").to_numpy()
other = ~ph

# Background trend: PH2 identity percentile as a function of cell-cycle percentile.
edges = np.linspace(0, 1, 21)
bin_id = np.clip(np.digitize(cycle_pct, edges) - 1, 0, 19)
centers = (edges[:-1] + edges[1:]) / 2
trend = np.full(20, np.nan)
for b in range(20):
    m = other & (bin_id == b)
    if m.sum() >= 40:
        trend[b] = np.median(prog_pct[m])
valid = np.isfinite(trend)
trend = np.interp(np.arange(20), np.flatnonzero(valid), trend[valid])
pred = np.interp(cycle_pct, centers, trend)
resid = prog_pct - pred
resid_z = (resid - resid.mean()) / resid.std(ddof=1)

# Plot
fig = plt.figure(figsize=(7.4, 4.5))
gs = fig.add_gridspec(1, 2, width_ratios=[1.55, 0.85], wspace=0.30, left=0.08, right=0.97, bottom=0.15, top=0.88)
ax = fig.add_subplot(gs[0, 0])
ax.hexbin(cycle_pct[other], prog_pct[other], gridsize=38, cmap="Greys", mincnt=1, linewidths=0)
ax.scatter(cycle_pct[ph], prog_pct[ph], s=4, alpha=0.22, color="#B2182B", edgecolors="none", rasterized=True)
# Background trend
ax.plot(centers, trend, color="#444444", lw=1.1)
# PH binned trend
ph_trend_x, ph_trend_y = [], []
for b in range(10):
    m = ph & (cycle_pct >= np.linspace(0, 1, 11)[b]) & (cycle_pct <= np.linspace(0, 1, 11)[b+1])
    if m.sum() >= 20:
        ph_trend_x.append(cycle_pct[m].mean()); ph_trend_y.append(prog_pct[m].mean())
ax.plot(ph_trend_x, ph_trend_y, color="#B2182B", lw=1.5)
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.set_xticks(np.linspace(0, 1, 6)); ax.set_yticks(np.linspace(0, 1, 6))
ax.set_xlabel("Cell-cycle program percentile")
ax.set_ylabel("PH identity program percentile")
ax.set_title("PH identity across cell-cycle states", pad=6)
ax.tick_params(length=3, width=0.6, pad=2)
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.text(0.03, 0.95, f"PH n={ph.sum():,}", transform=ax.transAxes,
        color="#B2182B", fontsize=7, weight="bold", va="top")
ax.text(0.03, 0.89, f"Other Hs cells n={other.sum():,}", transform=ax.transAxes,
        color="#444444", fontsize=6.5, va="top")

# Residual violin
ax2 = fig.add_subplot(gs[0, 1])
bg_vals = resid_z[other]
ph_vals = resid_z[ph]
vp1 = ax2.violinplot([bg_vals], positions=[0], widths=0.52, bw_method=0.35, points=220, showextrema=False)
vp2 = ax2.violinplot([ph_vals], positions=[1], widths=0.52, bw_method=0.35, points=220, showextrema=False)
for body in vp1["bodies"]:
    body.set_facecolor("#BDBDBD"); body.set_edgecolor("none"); body.set_alpha(0.72)
for body in vp2["bodies"]:
    body.set_facecolor("#B2182B"); body.set_edgecolor("none"); body.set_alpha(0.76)
ax2.plot([0, 1], [np.median(bg_vals), np.median(ph_vals)], color="#333333", lw=0.8, zorder=3)
ax2.scatter([0, 1], [np.median(bg_vals), np.median(ph_vals)], s=12, color=["#777777", "#B2182B"], zorder=4)
ax2.axhline(0, color="#777777", lw=0.7, ls="--")
ax2.set_ylim(-2.6, 2.6)
ax2.set_xticks([0, 1]); ax2.set_xticklabels(["Other", "PH"])
ax2.set_ylabel("Adjusted PH identity (z)")
ax2.set_title("Identity beyond cell cycle", pad=6)
ax2.text(1, 0.93, f"z={ph_vals.mean():.2f}", transform=ax2.get_xaxis_transform(), ha="center", va="top", color="#B2182B", fontsize=7, weight="bold")
ax2.tick_params(length=3, width=0.6)
for s in ("top", "right"): ax2.spines[s].set_visible(False)

fig.savefig(OUT / "Fig2D_replacement_preview.png", dpi=300, bbox_inches="tight")
fig.savefig(OUT / "Fig2D_replacement_preview.svg", format="svg", bbox_inches="tight")
plt.close(fig)
print(OUT / "Fig2D_replacement_preview.png")
