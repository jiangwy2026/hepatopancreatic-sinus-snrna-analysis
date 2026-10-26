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
t = pd.read_csv(out + r"\composition_pct_by_group.csv", index_col=0)
pct = t[["blank","treated"]].reindex(ORD)
W, THRESH, MINGAP = 0.44, 5.0, 3.8
X0, X1 = 0.0, 0.62
XMIN, XMAX = -0.70, 1.08
# 左侧引线终点（轴内），右侧引线终点（轴外页边）
XL_L, XL_R = -0.42, 0.93
BARS = [("blank", X0, "Control"), ("treated", X1, "6 h post-bleeding")]
fig = plt.figure(figsize=(5.0, 3.7))
ax = fig.add_axes([0.155, 0.13, 0.54, 0.775])
for g, x, lab in BARS:
    bot = 0.0
    for p in ORD:
        h = float(pct.loc[p, g])
        ax.bar(x, h, bottom=bot, width=W, color=COL[p], edgecolor="white", linewidth=0.4)
        if h >= THRESH:
            ax.text(x, bot + h/2, "%s\n%.1f%%" % (p, h), ha="center", va="center", fontsize=5.6,
                    weight="bold", color="black", linespacing=1.15)
        bot += h
for gi, (g, x, lab) in enumerate(BARS):
    side = -1 if gi == 0 else 1
    xline = XL_L if side < 0 else XL_R
    segs, bot = [], 0.0
    for p in ORD:
        h = float(pct.loc[p, g]); cy = bot + h/2; bot += h
        if h < THRESH: segs.append([p, h, cy])
    segs.sort(key=lambda s: -s[2])
    ys = [s[2] for s in segs]
    for i in range(1, len(ys)):
        if ys[i-1] - ys[i] < MINGAP: ys[i] = ys[i-1] - MINGAP
    for (p, h, cy), ly in zip(segs, ys):
        xedge = x + side*(W/2 + 0.015)
        xt = xline + side*0.03
        ax.plot([xedge, xline, xt], [cy, ly, ly], color=COL[p], lw=0.7,
                clip_on=False, solid_joinstyle="miter")
        ax.text(xt + side*0.02, ly, "%s %.1f%%" % (p, h), ha="left" if side > 0 else "right",
                va="center", fontsize=5.4, color=COL[p], weight="bold", clip_on=False)
ax.set_xlim(XMIN, XMAX); ax.set_ylim(0, 100)
ax.set_ylabel("Cells (%)", labelpad=2)
ax.set_xticks([X0, X1]); ax.set_xticklabels([b[2] for b in BARS], fontsize=6.8)
ax.tick_params(axis="x", length=2.5, width=0.6, pad=2)
ax.tick_params(axis="y", labelsize=6, length=2.2, width=0.6)
for sp in ("top","right"): ax.spines[sp].set_visible(False)
ax.spines["bottom"].set_visible(True)
ax.set_title("Cell composition in the hepatic sinus", fontsize=8.5, pad=7)
bp = os.path.join(figdir, "Fig3_composition_stacked_bar")
for _e, _kw in (("_preview.png", dict(dpi=300)), (".pdf", {}), (".svg", {}),
               (".tiff", dict(dpi=1200, pil_kwargs={"compression": "tiff_lzw"}))):
    fig.savefig(bp + _e, bbox_inches=None, **_kw)
    print("  %-52s %6.2f MB" % (os.path.basename(bp + _e), os.path.getsize(bp + _e) / 1e6))
plt.close()
print("预览已更新:", bp + "_preview.png")

