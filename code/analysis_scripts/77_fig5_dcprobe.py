import sys, io, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd, scanpy as sc
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 250); sc.settings.verbosity = 0
out = r"D:\桌面\single cell\analysis_output"; figdir = out + r"\figures_pub"
H = sc.read_h5ad(out + r"\adata_hemocyte_fig5.h5ad")
dc = H.obsm["X_diffmap"]; pt = H.obs["dpt_pseudotime"].values.astype(float)
print("扩散分量离散度:", [round(float(dc[:,i].std()),4) for i in range(6)])
combos = [(1,2), (0,1), (2,3), (0,2)]
fig, axes = plt.subplots(2, len(combos), figsize=(4.0*len(combos), 7.4))
for j,(a,b) in enumerate(combos):
    for i,(ax, mode) in enumerate(zip(axes[:,j], ["pt","ct"])):
        if mode == "pt":
            v = pt.copy(); v[~np.isfinite(v)] = np.nan
            s = ax.scatter(dc[:,a], dc[:,b], c=v, s=1.4, cmap="Blues", linewidths=0, rasterized=True)
            plt.colorbar(s, ax=ax, fraction=0.046, pad=0.02)
        else:
            for c,col in (("PH","#B2182B"),("GH","#E08214"),("MH1","#FDB863"),("MH2","#8C510A")):
                m = (H.obs["ct"]==c).values
                ax.scatter(dc[m,a], dc[m,b], s=1.4, c=col, linewidths=0, rasterized=True, label=c)
            ax.legend(fontsize=5, markerscale=4, frameon=False)
        ax.set_title("DC%d vs DC%d (%s)" % (a,b,mode), fontsize=7)
        ax.set_xticks([]); ax.set_yticks([])
plt.tight_layout()
fig.savefig(os.path.join(figdir, "_fig5_dc_options.png"), dpi=150)
plt.close()
print("已出对比图: _fig5_dc_options.png")
