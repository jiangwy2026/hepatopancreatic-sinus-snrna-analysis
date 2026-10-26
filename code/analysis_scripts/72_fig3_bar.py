import sys, io, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
out = r"D:\桌面\single cell\analysis_output"; figdir = out + r"\figures_pub"
ORD = ["PH","GH","MH1","MH2","Epi","F","R","B","NEU"]
COL = {"PH":"#B2182B","GH":"#E08214","MH1":"#FDB863","MH2":"#8C510A","Epi":"#7B3294",
       "F":"#006D2C","R":"#41AB5D","B":"#78C679","NEU":"#2166AC"}
cnt = pd.crosstab(pd.read_csv(out + r"\composition_pct_by_group.csv", index_col=0).index, 1)  # placeholder
t = pd.read_csv(out + r"\composition_pct_by_group.csv", index_col=0)
pct = t[["blank","treated"]].reindex(ORD)
fc  = t["倍数"].reindex(ORD)
gtot = pd.Series({"blank": 5349, "treated": 4309})
fig = plt.figure(figsize=(7.4, 3.5))
gs = fig.add_gridspec(1, 3, width_ratios=[0.62, 1.0, 1.15], wspace=0.45, left=0.075, right=0.985, top=0.90, bottom=0.13)
# ---- (a) 堆叠柱状图 ----
axa = fig.add_subplot(gs[0])
groups = ["blank","treated"]
labels = ["Control", "6 h post-bleeding"]
bot = np.zeros(2)
for g in ORD:
    v = pct.loc[g, groups].values.astype(float)
    axa.bar([0,1], v, bottom=bot, color=COL[g], width=0.62, edgecolor="white", linewidth=0.35, label=g)
    for xi, (b, h) in enumerate(zip(bot, v)):
        if h >= 5:
            axa.text(xi, b + h/2, "%s\n%.1f%%" % (g, h), ha="center", va="center", fontsize=5.2,
                     color="black", weight="bold")
    bot += v
axa.set_xticks([0,1]); axa.set_xticklabels(labels, fontsize=6.5)
axa.set_ylim(0, 100); axa.set_ylabel("Cells (%)", labelpad=2)
axa.set_title("Composition", fontsize=8, pad=4)
axa.tick_params(axis="both", labelsize=6, length=2.2, width=0.6)
for sp in ("top","right"): axa.spines[sp].set_visible(False)
# ---- (b) 各亚群占比对比 ----
axb = fig.add_subplot(gs[1])
y = np.arange(len(ORD))[::-1]; h = 0.36
axb.barh(y + h/2, pct["blank"].values,   height=h, color="#BFBFBF", label="Control", edgecolor="none")
axb.barh(y - h/2, pct["treated"].values, height=h, color="#C0392B", label="6 h post-bleeding", edgecolor="none")
axb.set_yticks(y); axb.set_yticklabels(ORD, fontsize=6.5)
axb.set_xlabel("Cells per group (%)", labelpad=2); axb.set_title("Population shift", fontsize=8, pad=4)
axb.tick_params(axis="both", labelsize=6, length=2.2, width=0.6)
axb.legend(fontsize=5.8, loc="lower right", frameon=False, handlelength=1.0, handletextpad=0.4)
for sp in ("top","right"): axb.spines[sp].set_visible(False)
# ---- (c) 倍数变化 ----
axc = fig.add_subplot(gs[2])
vals = np.log2(fc.values)
cols = ["#C0392B" if v > 0 else "#4C72B0" for v in vals]
axc.barh(y, vals, height=0.62, color=cols, edgecolor="none")
axc.axvline(0, color="black", lw=0.6)
axc.set_yticks(y); axc.set_yticklabels([])
axc.set_xlabel("log2 fold change (6 h / control)", labelpad=2); axc.set_title("Enrichment after bleeding", fontsize=8, pad=4)
axc.tick_params(axis="both", labelsize=6, length=2.2, width=0.6)
for sp in ("top","right"): axc.spines[sp].set_visible(False)
for yi, (g, v, f) in zip(y, zip(ORD, vals, fc.values)):
    axc.text(v + (0.12 if v >= 0 else -0.12), yi, "%.2f×" % f, va="center",
             ha="left" if v >= 0 else "right", fontsize=5.3, color="black")
axc.set_xlim(min(vals)*1.45, max(vals)*1.45)
fig.text(0.012, 0.965, "a", fontsize=9, weight="bold"); fig.text(0.245, 0.965, "b", fontsize=9, weight="bold"); fig.text(0.60, 0.965, "c", fontsize=9, weight="bold")
bp = os.path.join(figdir, "Fig3_composition_stacked_bar")
for e, kw in (("_preview.png", dict(dpi=300)), (".pdf", {}), (".svg", {}), (".tiff", dict(dpi=1200, pil_kwargs={"compression":"tiff_lzw"}))):
    fig.savefig(bp + e, bbox_inches="tight", pad_inches=0.04, **kw)
plt.close()
for e in ("_preview.png",".pdf",".svg",".tiff"): print("  %-38s %6.2f MB" % (os.path.basename(bp+e), os.path.getsize(bp+e)/1e6))
