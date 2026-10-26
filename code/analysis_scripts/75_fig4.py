import sys, io, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np, pandas as pd, scanpy as sc
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 250); sc.settings.verbosity = 0
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"; figdir = out + r"\figures_pub"
H = sc.read_h5ad(out + r"\adata_hemocyte_fig4.h5ad")
H.X = H.layers["counts"].copy(); sc.pp.normalize_total(H, target_sum=1e4); sc.pp.log1p(H)
sym2i = {}
for i, s in zip(H.var_names, H.var["sym"].astype(str)): sym2i.setdefault(s, i)
pt = H.obs["dpt_pseudotime"].values.astype(float)
finite = np.isfinite(pt)
print("拟时序: 有限 %d / %d (排除 %d 个 inf/NaN)" % (finite.sum(), len(pt), (~finite).sum()))
ptf = pt[finite]
nb = 20                                    # 等细胞数分箱（分位数），避免拟时序偏斜导致空箱
edges = np.quantile(ptf, np.linspace(0, 1, nb+1))
edges[-1] += 1e-9
centers = np.array([ptf[(ptf >= edges[i]) & (ptf < edges[i+1])].mean() for i in range(nb)])
binid = np.full(len(pt), -1)
binid[finite] = np.clip(np.digitize(ptf, edges)-1, 0, nb-1)
def binned(v):
    m = np.array([v[finite & (binid == b)].mean() if (finite & (binid == b)).sum() else np.nan for b in range(nb)])
    return pd.Series(m).rolling(3, center=True, min_periods=1).mean().values   # 轻度平滑
COL = {"PH":"#B2182B","GH":"#E08214","MH1":"#FDB863"}
um = H.obsm["X_umap"]
fig = plt.figure(figsize=(7.4, 6.3))
gs = fig.add_gridspec(2, 2, wspace=0.30, hspace=0.36, left=0.085, right=0.975, top=0.94, bottom=0.075)
# (a)
axa = fig.add_subplot(gs[0,0])
ptc = pt.copy(); ptc[~finite] = np.nan
vmax = np.nanpercentile(ptc, 99)
sc0 = axa.scatter(um[:,0], um[:,1], c=ptc, s=1.6, cmap='viridis', vmin=0, vmax=vmax, linewidths=0, rasterized=True)
iroot = int(np.nanargmin(pt))
axa.scatter(um[iroot,0], um[iroot,1], s=45, marker="*", c="white", edgecolors="black", linewidths=0.6, zorder=5)
axa.annotate("root (PH)", xy=(um[iroot,0], um[iroot,1]), xytext=(um[iroot,0]-5.2, um[iroot,1]+3.6),
             fontsize=6, arrowprops=dict(arrowstyle="->", lw=0.6, color="black"))
cb = fig.colorbar(sc0, ax=axa, fraction=0.045, pad=0.02); cb.set_label("Pseudotime", fontsize=6)
cb.ax.tick_params(labelsize=5.5); cb.outline.set_visible(False)
axa.set_title("Hemocyte trajectory", fontsize=8, pad=4)
# (b)
axb = fig.add_subplot(gs[0,1])
for c in ["PH","GH","MH1"]:
    m = (H.obs["cell"] == c).values
    axb.scatter(um[m,0], um[m,1], s=1.6, c=COL[c], linewidths=0, rasterized=True, label="%s (%d)" % (c, m.sum()))
axb.legend(fontsize=6, markerscale=4, loc="upper right", frameon=False, handletextpad=0.3, labelspacing=0.3)
axb.set_title("Hemocyte subtypes", fontsize=8, pad=4)
for ax in (axa, axb):
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP 1", labelpad=1); ax.set_ylabel("UMAP 2", labelpad=1)
    for sp in ("top","right","left","bottom"): ax.spines[sp].set_visible(False)
# (c)
axc = fig.add_subplot(gs[1,0])
prog = ["LOC113823945","ccnd2","Igfbp7","Grik2"]; matu = ["LOC113801825","PEN-3","Pxt","A2M"]
reds = ["#67000D","#A50F15","#CB181D","#EF3B2C"]; blues = ["#084594","#2171B5","#4292C6","#6BAED6"]
for g, c in zip(prog, reds):
    if g not in sym2i: continue
    v = np.asarray(H[:, sym2i[g]].X.todense()).ravel(); m = binned(v)
    m = (m-np.nanmin(m))/(np.nanmax(m)-np.nanmin(m)+1e-9)
    axc.plot(centers, m, color=c, lw=1.1, label=g)
for g, c in zip(matu, blues):
    if g not in sym2i: continue
    v = np.asarray(H[:, sym2i[g]].X.todense()).ravel(); m = binned(v)
    m = (m-np.nanmin(m))/(np.nanmax(m)-np.nanmin(m)+1e-9)
    axc.plot(centers, m, color=c, lw=1.1, label=g)
axc.set_xlabel("Pseudotime", labelpad=2); axc.set_ylabel("Relative expression", labelpad=2)
axc.set_ylim(-0.03, 1.15); axc.set_title("Marker dynamics", fontsize=8, pad=4)
axc.legend(fontsize=5.0, ncol=2, frameon=False, loc="lower right", handlelength=1.2, handletextpad=0.35, labelspacing=0.28, columnspacing=0.7)
axc.tick_params(labelsize=6, length=2.2, width=0.6)
for sp in ("top","right"): axc.spines[sp].set_visible(False)
# (d)
axd = fig.add_subplot(gs[1,1])
P = pd.read_csv(out + r"\marker_heatmap_selected.csv")
for grp, c, lab in (("PH","#B2182B","Progenitor (PH)"), ("GH","#2171B5","Mature effector (GH)")):
    gs_ = [g for g in P[P.grp == grp].gene if g in sym2i]
    v = np.asarray(H[:, [sym2i[g] for g in gs_]].X.mean(1)).ravel(); m = binned(v)
    m = (m-np.nanmin(m))/(np.nanmax(m)-np.nanmin(m)+1e-9)
    axd.plot(centers, m, color=c, lw=1.5, label=lab)
axd.set_xlabel("Pseudotime", labelpad=2); axd.set_ylabel("Programme score (scaled)", labelpad=2)
axd.set_ylim(-0.03, 1.15); axd.set_title("Programme switch", fontsize=8, pad=4)
axd.legend(fontsize=5.8, frameon=False, loc="center right", handlelength=1.4, handletextpad=0.4)
axd.tick_params(labelsize=6, length=2.2, width=0.6)
for sp in ("top","right"): axd.spines[sp].set_visible(False)
for ax, lt in ((axa,"a"),(axb,"b"),(axc,"c"),(axd,"d")):
    ax.text(-0.06, 1.06, lt, transform=ax.transAxes, fontsize=9, weight="bold", va="bottom")
bp = os.path.join(figdir, "Fig4_pseudotime")
fig.savefig(bp + "_preview.png", dpi=300, bbox_inches=None)
plt.close()
print("预览:", bp + "_preview.png")


