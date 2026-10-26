import sys, io, os, traceback
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.gridspec as gs
import numpy as np, pandas as pd, scanpy as sc
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style, save_pub
apply_pub_style(7)
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"; figdir = out + r"\figures_pub"; os.makedirs(figdir, exist_ok=True)
PAL = ["#4C72B0","#DD8452","#55A868","#C44E52","#8172B3","#937860","#DA8BC3","#8C8C8C","#CCB974","#64B5CD"]
def save(fig, name, **kw):
    p = save_pub(fig, os.path.join(figdir, name), **kw); plt.close(fig)
    print("  saved:", [os.path.basename(x) for x in p])

# ---------- Fig 1: 全局概览 ----------
try:
    G = sc.read_h5ad(out + r"\adata_hvig_processed.h5ad") if False else sc.read_h5ad(out + r"\adata_hvg_processed.h5ad")
    cm = pd.read_csv(out + r"\cellmeta_raw.csv")
    fig = plt.figure(figsize=(7.2, 2.5))
    g = gs.GridSpec(1, 3, width_ratios=[1.15, 1.15, 0.95], wspace=0.35, figure=fig)
    ax0 = fig.add_subplot(g[0]); ax1 = fig.add_subplot(g[1]); ax2 = fig.add_subplot(g[2])
    X = G.obsm["X_umap"]; cl = G.obs["vendor_cluster"].astype(str).values
    order = sorted(set(cl), key=int)
    for i, c in enumerate(order):
        m = cl == c
        ax0.scatter(X[m,0], X[m,1], s=1.6, c=PAL[i % len(PAL)], linewidths=0, rasterized=True, label=c)
        ax0.text(np.median(X[m,0]), np.median(X[m,1]), c, fontsize=6, ha="center", va="center", weight="bold")
    ax0.set_title("Clusters (n = 10)", fontsize=8)
    samp = G.obs["sample"].astype(str).values if "sample" in G.obs else np.array(["?"]*len(cl))
    for i, s in enumerate(["blank","treated"]):
        m = samp == s
        ax1.scatter(X[m,0], X[m,1], s=1.6, c=["#9aa4b0","#C44E52"][i], linewidths=0, rasterized=True, label=s)
    ax1.set_title("Sample", fontsize=8); ax1.legend(markerscale=6, loc="best")
    for a in (ax0, ax1):
        a.set_xticks([]); a.set_yticks([]); a.set_xlabel("UMAP1"); a.set_ylabel("UMAP2")
        for sp in ("top","right","left","bottom"): a.spines[sp].set_visible(False)
    t = pd.crosstab(cm.Cluster, cm.Groups)
    t = t.div(t.sum(0), axis=1) * 100
    y = np.arange(len(t)); h = 0.38
    ax2.barh(y + h/2, t["blank"], height=h, color="#9aa4b0", label="blank")
    ax2.barh(y - h/2, t["treated"], height=h, color="#C44E52", label="treated")
    ax2.set_yticks(y); ax2.set_yticklabels([f"C{i}" for i in t.index]); ax2.invert_yaxis()
    ax2.set_xlabel("Cells (%)"); ax2.set_title("Composition", fontsize=8); ax2.legend(loc="lower right")
    for i, c in enumerate(t.index):
        fc = t.loc[c, "treated"] / t.loc[c, "blank"]
        if fc >= 1.8 or fc <= 0.5:
            ax2.text(max(t.loc[c,"blank"], t.loc[c,"treated"]) + 2, i, f"{fc:.1f}×", fontsize=6, va="center")
    fig.suptitle("Hepatopancreatic sinus scRNA-seq: 10 clusters, control vs 6 h post-bleeding", fontsize=8, y=1.04)
    save(fig, "Fig1_UMAP_overview")
except Exception: traceback.print_exc()

# ---------- Fig 2: 血细胞轨迹 ----------
try:
    H = sc.read_h5ad(out + r"\adata_hemocyte_lineage.h5ad")
    sub = H.obs["sub"].astype(str).values; ph = H.obs["phase"].astype(str).values
    pt = H.obs["dpt_pseudotime"].values
    fate = H.obs["fate_to_1"].values if "fate_to_1" in H.obs else np.full(len(sub), np.nan)
    X = H.obsm["X_umap"]
    fig = plt.figure(figsize=(7.2, 2.1))
    g = gs.GridSpec(1, 4, wspace=0.25, figure=fig)
    specs = [("Hemocyte subtypes", None), ("Cell-cycle phase", None), ("Pseudotime (root = cycling)", pt), ("Fate prob. to granular", fate)]
    for k, (title, val) in enumerate(specs):
        ax = fig.add_subplot(g[k])
        if val is None and k == 0:
            for i, s in enumerate(sorted(set(sub))):
                m = sub == s
                ax.scatter(X[m,0], X[m,1], s=1.4, c=PAL[i % len(PAL)], linewidths=0, rasterized=True)
                ax.text(np.median(X[m,0]), np.median(X[m,1]), s, fontsize=6, ha="center", va="center", weight="bold")
        elif val is None and k == 1:
            for i, s in enumerate(["G1","S","G2M"]):
                m = ph == s
                ax.scatter(X[m,0], X[m,1], s=1.4, c=["#cccccc","#8172B3","#C44E52"][i], linewidths=0, rasterized=True, label=s)
            ax.legend(markerscale=6, loc="best")
        else:
            o = np.argsort(np.nan_to_num(val, nan=-1))
            ax.scatter(X[o,0], X[o,1], s=1.4, c=val[o], cmap="viridis", linewidths=0, rasterized=True)
        ax.set_title(title, fontsize=7.5); ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP1")
        for sp in ("top","right","left","bottom"): ax.spines[sp].set_visible(False)
    fig.suptitle("Cycling hemocytes (subtype 2) give rise to mature effector hemocytes", fontsize=8, y=1.06)
    save(fig, "Fig2_hemocyte_trajectory")
except Exception: traceback.print_exc()

# ---------- Fig 3: PAGA 全局 ----------
try:
    Gg = sc.read_h5ad(out + r"\adata_global_integrated.h5ad")
    fig, ax = plt.subplots(figsize=(3.4, 3.2))
    sc.pl.paga(Gg, color="grp", ax=ax, show=False, node_size_scale=2.2, edge_width_scale=3.5, fontsize=6, legend_loc=None)
    ax.set_title("Global PAGA topology", fontsize=8)
    save(fig, "Fig3_PAGA_global")
except Exception: traceback.print_exc()

# ---------- Fig 4: marker dotplot ----------
try:
    Hh = sc.read_h5ad(out + r"\adata_hemocyte_harmony.h5ad")
    A = sc.read_h5ad(out + r"\adata_raw.h5ad"); A = A[list(Hh.obs_names)].copy()
    A.X = A.layers["counts"].copy(); sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
    A.obs["sub"] = Hh.obs["sub_0.3"].astype(str).values
    A.var_names = [str(x) for x in A.var["symbol"].values]; A.var_names_make_unique()
    panel = ["A2M","PPAF2","Pxt","FCN2","pigs","ccnd2","LOC113823945","SSPO","Grik2","slo",
             "THBS4","Act5C","eef1a","RpLP1","Apod","LOC113824714","Crabp2","LOC113828299"]
    panel = [g for g in panel if g in set(A.var_names)]
    ax = sc.pl.dotplot(A, panel, groupby="sub", show=False, return_fig=False)
    fig = plt.gcf(); fig.set_size_inches(6.6, 2.0)
    fig.axes[0].set_title("Hemocyte subtype markers", fontsize=8)
    save(fig, "Fig4_marker_dotplot")
except Exception: traceback.print_exc()

print("\n=== figures_pub 目录 ===")
for f in sorted(os.listdir(figdir)):
    print("  %-34s %8.2f MB" % (f, os.path.getsize(os.path.join(figdir, f))/1e6))
