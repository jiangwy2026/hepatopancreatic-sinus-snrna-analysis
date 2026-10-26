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
fig, ax = plt.subplots(figsize=(4.6, 4.2))
for ab in order:
    m = (A.obs["ct"] == ab).values
    ax.scatter(um[m,0], um[m,1], s=2.0, c=COL[ab], linewidths=0, rasterized=True,
               label="%s  %s (%d)" % (ab, FULL[ab], int(m.sum())))
for ab in order:
    m = (A.obs["ct"] == ab).values
    cx, cy = np.median(um[m,0]), np.median(um[m,1])
    if ab == "NEU": cx += 1.2
    ax.text(cx, cy, ab, fontsize=7, weight="bold", ha="center", va="center",
            path_effects=[pe.withStroke(linewidth=2.2, foreground="white")])
ax.set_aspect("equal")
ax.set_xticks(np.arange(-5, 16, 5)); ax.set_yticks(np.arange(0, 18, 2))
ax.set_xlabel("UMAP 1", labelpad=2); ax.set_ylabel("UMAP 2", labelpad=2)
ax.tick_params(axis="both", labelsize=6, length=2.2, width=0.6)
for sp in ("top","right"): ax.spines[sp].set_visible(False)
ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=6.0, handletextpad=0.35,
          labelspacing=0.42, borderaxespad=0.0, markerscale=4, frameon=False)
ax.set_title("Hepatic sinus cell populations (n = 9,658)", fontsize=8.5, pad=6)
plt.tight_layout()
bp = os.path.join(figdir, "Fig1_UMAP_celltypes")
fig.savefig(bp + "_preview.png", dpi=600, bbox_inches="tight", pad_inches=0.03)
fig.savefig(bp + ".pdf", bbox_inches="tight", pad_inches=0.03)
fig.savefig(bp + ".svg", bbox_inches="tight", pad_inches=0.03)
fig.savefig(bp + ".tiff", dpi=1200, bbox_inches="tight", pad_inches=0.03,
            pil_kwargs={"compression": "tiff_lzw"})
plt.close()
for e in ("_preview.png", ".pdf", ".svg", ".tiff"):
    f = bp + e
    print("  %-30s %7.2f MB" % (os.path.basename(f), os.path.getsize(f)/1e6))
