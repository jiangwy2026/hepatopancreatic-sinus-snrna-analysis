import sys, io, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np, pandas as pd, scanpy as sc
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"; figdir = out + r"\figures_pub"
A = sc.read_h5ad(out + r"\adata_soupx_hvg.h5ad")
base = {"2":"PH","1":"GH","0":"MH1","3":"MH2","4":"Epi","5":"F","6":"R","7":"B","8":"NE","9":"Neu"}
A.obs["ct"] = A.obs.vendor_cluster.astype(str).map(base).replace({"NE":"NEU","Neu":"NEU"})
FULL = {"PH":"Prohemocyte","GH":"Granular hemocyte","MH1":"Monocytic hemocyte 1","MH2":"Macrophage-like phagocyte",
        "Epi":"Epidermal / cuticle cell","F":"Fibrillar cell","R":"Resorptive cell","B":"Blister-like cell",
        "NEU":"Neural / neuroendocrine cell"}
COL = {"PH":"#B2182B","GH":"#E08214","MH1":"#FDB863","MH2":"#8C510A","Epi":"#7B3294",
       "F":"#006D2C","R":"#41AB5D","B":"#78C679","NEU":"#2166AC"}
order = ["PH","GH","MH1","MH2","Epi","F","R","B","NEU"]
um = A.obsm["X_umap"]
def draw(ax, style):
    for ab in order:
        m = (A.obs["ct"] == ab).values
        ax.scatter(um[m,0], um[m,1], s=2.0, c=COL[ab], linewidths=0, rasterized=True,
                   label="%s  %s (%d)" % (ab, FULL[ab], int(m.sum())))
    for ab in order:
        m = (A.obs["ct"] == ab).values
        cx, cy = np.median(um[m,0]), np.median(um[m,1])
        if ab == "NEU": cx += 1.2
        kw = dict(fontsize=7, weight="bold", ha="center", va="center")
        if style == "outline":
            kw["path_effects"] = [pe.withStroke(linewidth=2.2, foreground="white")]
        ax.text(cx, cy, ab, **kw)
    ax.set_aspect("equal")
    ax.set_xticks(np.arange(-5, 16, 5)); ax.set_yticks(np.arange(0, 18, 2))
    ax.set_xlabel("UMAP 1", labelpad=2); ax.set_ylabel("UMAP 2", labelpad=2)
    ax.tick_params(axis="both", labelsize=6, length=2.2, width=0.6)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)
fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.2))
draw(axes[0], "plain");   axes[0].set_title("A  Plain text (no box, no outline)", fontsize=8)
draw(axes[1], "outline"); axes[1].set_title("B  Text with white outline (no box)", fontsize=8)
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, loc="center right", bbox_to_anchor=(1.005, 0.5), fontsize=6.0,
           handletextpad=0.35, labelspacing=0.42, markerscale=4, frameon=False)
plt.tight_layout(rect=(0, 0, 0.80, 1))
p = os.path.join(figdir, "Fig1_label_style_options.png")
fig.savefig(p, dpi=300, bbox_inches="tight", pad_inches=0.04); plt.close()
print("ok:", p)
