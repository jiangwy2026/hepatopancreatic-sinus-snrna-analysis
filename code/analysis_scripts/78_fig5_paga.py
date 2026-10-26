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
pt = H.obs["dpt_pseudotime"].values.astype(float)
sc.tl.paga(H, groups='ct')
sc.pl.paga(H, show=False, plot=False)
ct = pd.DataFrame(H.uns["paga"]["connectivities"].todense(), index=H.obs["ct"].cat.categories, columns=H.obs["ct"].cat.categories)
print("=== PAGA 连通性 ==="); print(ct.round(3).to_string())
# PAGA 引导的布局
try:
    sc.tl.draw_graph(H, init_pos="paga", random_state=0)
    print("draw_graph(paga init) ok")
except Exception as e:
    print("draw_graph 失败:", str(e)[:100])
sc.tl.umap(H, init_pos="paga", random_state=0)
fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
for ax, key, ttl in ((axes[0], "X_draw_graph_fa", "draw_graph (PAGA init)"), (axes[1], "X_umap", "UMAP (PAGA init)")):
    if key not in H.obsm: continue
    X = H.obsm[key]
    v = pt.copy(); v[~np.isfinite(v)] = np.nan
    s = ax.scatter(X[:,0], X[:,1], c=v, s=1.4, cmap="Blues", linewidths=0, rasterized=True)
    plt.colorbar(s, ax=ax, fraction=0.046, pad=0.02)
    ax.set_title(ttl + " — pseudotime", fontsize=7)
ax3 = axes[2]
if "X_umap" in H.obsm:
    X = H.obsm["X_umap"]
    for c, col in (("PH","#B2182B"),("GH","#E08214"),("MH1","#FDB863"),("MH2","#8C510A")):
        m = (H.obs["ct"]==c).values
        ax3.scatter(X[m,0], X[m,1], s=1.4, c=col, linewidths=0, rasterized=True, label=c)
    ax3.legend(fontsize=5.5, markerscale=4, frameon=False)
    ax3.set_title("UMAP (PAGA init) — subtypes", fontsize=7)
for ax in axes: ax.set_xticks([]); ax.set_yticks([])
plt.tight_layout(); fig.savefig(os.path.join(figdir, "_fig5_paga_options.png"), dpi=150); plt.close()
H.write(out + r"\adata_hemocyte_fig5b.h5ad")
print("已出对比: _fig5_paga_options.png")

