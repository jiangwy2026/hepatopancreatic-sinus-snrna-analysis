import sys, io, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Ellipse
import numpy as np, pandas as pd, scanpy as sc
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 250); sc.settings.verbosity = 0
out = r"D:\桌面\single cell\analysis_output"; figdir = out + r"\figures_pub"
H = sc.read_h5ad(out + r"\adata_hemocyte_fig5.h5ad")     # 含 X_pca_harmony + dpt_pseudotime
sc.pp.neighbors(H, use_rep="X_pca_harmony", n_neighbors=15, n_pcs=30)
sc.tl.umap(H, random_state=0, init_pos="spectral")
H.write(out + r"\adata_hemocyte_fig5c.h5ad")
X = H.obsm["X_umap"]; pt = H.obs["dpt_pseudotime"].values.astype(float)
ct = H.obs["ct"].astype(str).values
COL = {"PH":"#B2182B","GH":"#2C7FB8","MH1":"#F0A030","MH2":"#8C510A"}
LIN = [("Prohemocyte", ["PH"], "#B2182B"), ("Granulocyte lineage", ["GH"], "#2C7FB8"),
       ("Monocytic hemocyte lineage", ["MH1","MH2"], "#8C510A")]
fig = plt.figure(figsize=(9.6, 4.6))
axa = fig.add_axes([0.045, 0.055, 0.355, 0.815])
axb = fig.add_axes([0.515, 0.055, 0.355, 0.815])
v = pt.copy(); v[~np.isfinite(v)] = np.nan
s = axa.scatter(X[:,0], X[:,1], c=v, s=1.7, cmap="Blues_r", vmin=0, vmax=np.nanpercentile(v,99), linewidths=0, rasterized=True)
iroot = int(np.nanargmin(pt))
axa.scatter(X[iroot,0], X[iroot,1], s=60, marker="*", c="#E31A1C", edgecolors="white", linewidths=0.7, zorder=8)
axa.annotate("root (PH)", xy=(X[iroot,0], X[iroot,1]),
             xytext=(X[iroot,0]-6, X[iroot,1]+4), fontsize=6.5, weight="bold",
             path_effects=[pe.withStroke(linewidth=2.2, foreground="white")],
             arrowprops=dict(arrowstyle="->", lw=0.7, color="black"))
cax = fig.add_axes([0.075, 0.935, 0.22, 0.028])
cb = fig.colorbar(s, cax=cax, orientation="horizontal"); cb.ax.tick_params(labelsize=5.5, length=1.5); cb.outline.set_visible(False)
cb.set_label("Pseudotime", fontsize=6.2, labelpad=2)

for c in ["PH","GH","MH1","MH2"]:
    m = ct == c
    axb.scatter(X[m,0], X[m,1], s=1.7, c=COL[c], linewidths=0, rasterized=True, label="%s (%d)" % (c, m.sum()))
for lab, groups, col in LIN:
    m = np.isin(ct, groups); P = X[m]
    mu = P.mean(0); C = np.cov(P.T); w, vec = np.linalg.eigh(C)
    ang = np.degrees(np.arctan2(vec[1,-1], vec[0,-1]))
    axb.add_patch(Ellipse(mu, 2*1.75*np.sqrt(w[-1]), 2*1.75*np.sqrt(w[0]), angle=ang,
                          fill=False, ls=(0,(5,2.5)), lw=1.0, edgecolor=col, zorder=6))
    dy = {"Prohemocyte": 0.95, "Granulocyte lineage": 0.95, "Monocytic hemocyte lineage": -1.05}[lab]
    axb.text(mu[0], mu[1] + dy*1.75*np.sqrt(w[-1])*0.75, lab, fontsize=6.0, color=col, ha="center",
             va="bottom" if dy > 0 else "top", zorder=7, weight="bold",
             path_effects=[pe.withStroke(linewidth=2.6, foreground="white")])
axb.legend(fontsize=6, markerscale=4, frameon=False, loc="lower left", handletextpad=0.3, labelspacing=0.3)
axb.set_title("Subtypes and lineages", fontsize=8.5, pad=6)
for ax in (axa, axb):
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP 1", labelpad=1); ax.set_ylabel("UMAP 2", labelpad=1)
    for sp in ("top","right","left","bottom"): ax.spines[sp].set_visible(False)
fig.text(0.012, 0.955, "D.", fontsize=10, weight="bold")
fig.text(0.482, 0.955, "E.", fontsize=10, weight="bold")
bp = os.path.join(figdir, "Fig5_trajectory_lineages")
fig.savefig(bp + "_preview.png", dpi=300, bbox_inches=None)
plt.close()
print("预览:", bp + "_preview.png")

