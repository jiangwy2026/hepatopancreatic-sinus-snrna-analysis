import sys, io, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import numpy as np, pandas as pd, scanpy as sc
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 250); sc.settings.verbosity = 0
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"; figdir = out + r"\figures_pub"
A = sc.read_h5ad(out + r"\adata_soupx.h5ad")
base = {"2":"PH","1":"GH","0":"MH1","3":"MH2","4":"Epi","5":"F","6":"R","7":"B","8":"NE","9":"Neu"}
A.obs["ct"] = A.obs.vendor_cluster.astype(str).map(base).replace({"NE":"NEU","Neu":"NEU"})
ORD = ["PH","GH","MH1","MH2","Epi","F","R","B","NEU"]
A.obs["ct"] = pd.Categorical(A.obs["ct"], categories=ORD)
A.X = A.layers["counts"].copy(); sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
COL = {"PH":"#B2182B","GH":"#E08214","MH1":"#FDB863","MH2":"#8C510A","Epi":"#7B3294",
       "F":"#006D2C","R":"#41AB5D","B":"#78C679","NEU":"#2166AC"}
mt = {"ND1","ND2","ND3","ND4","ND4L","ND5","COX1","COX2","COX3","ATP6","ATP8","CYTB"}
symv = A.var["sym"].astype(str).values
sc.tl.rank_genes_groups(A, "ct", method="wilcoxon", pts=True)
mk = pd.concat([sc.get.rank_genes_groups_df(A, group=g).assign(grp=g) for g in ORD], ignore_index=True)
mk["gene"] = mk["names"].map(dict(zip(A.var_names, symv))).fillna(mk["names"]).astype(str)
mk = mk[~mk.gene.isin(mt)].copy()
mk["named"] = ~mk["gene"].str.startswith(("LOC","gene-LOC"))
mk["det_diff"] = mk.pct_nz_group - mk.pct_nz_reference
pool = mk[(mk.pct_nz_group > 0.15) & (mk.det_diff > 0.05) & (mk.logfoldchanges > 1.0)].copy()
pool["rank_score"] = pool.logfoldchanges * pool.det_diff
pool = pool.sort_values(["named","rank_score"], ascending=[False, False]).drop_duplicates("gene")
picked = []
for g in ORD:
    for _, r in pool[pool.grp == g].head(5).iterrows():
        picked.append({"grp": g, "gene": r.gene})
P = pd.DataFrame(picked)
print("每群:", dict(P.grp.value_counts().reindex(ORD)))
print(P.groupby("grp", sort=False).gene.apply(lambda s: ", ".join(s)).to_string())
sym2i = {}
for i, s in zip(A.var_names, symv): sym2i.setdefault(s, i)
Z = np.zeros((len(P), len(ORD)))
for r_i, g in enumerate(P.gene):
    j = sym2i[g]
    for c_i, gg in enumerate(ORD):
        m = (A.obs["ct"] == gg).values
        Z[r_i, c_i] = float(np.asarray(A[m][:, j].X.mean()).ravel()[0])
Zz = (Z - Z.mean(1, keepdims=True)) / (Z.std(1, keepdims=True) + 1e-9)
P.to_csv(out + r"\marker_heatmap_selected.csv", index=False)
pd.DataFrame(Zz, index=P.gene, columns=ORD).to_csv(out + r"\marker_heatmap_z.csv")
# ---------- 绘图：4 个轴（标签 / 色条 / 热图 / 色标）----------
n = len(P)
fig = plt.figure(figsize=(4.7, 8.4))
gs = fig.add_gridspec(1, 5, width_ratios=[1.15, 0.20, 0.09, 1.0, 0.05], wspace=0.05,
                      left=0.03, right=0.93, top=0.94, bottom=0.04)
axl = fig.add_subplot(gs[0]); axg = fig.add_subplot(gs[1]); axs = fig.add_subplot(gs[2]); axh = fig.add_subplot(gs[3]); axc = fig.add_subplot(gs[4])
gidx = np.array([ORD.index(g) for g in P.grp]).reshape(-1, 1)
axs.imshow(gidx, aspect="auto", cmap=ListedColormap([COL[g] for g in ORD]),
           norm=BoundaryNorm(np.arange(-0.5, len(ORD)+0.5), len(ORD)), extent=(0,1,n-0.5,-0.5))
axs.set_xticks([]); axs.set_yticks([])
for sp in axs.spines.values(): sp.set_visible(False)
im = axh.imshow(Zz, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
axh.set_xticks(range(len(ORD))); axh.set_xticklabels(ORD, rotation=45, ha="left", fontsize=6.5)
axh.xaxis.set_ticks_position("top"); axh.set_yticks([]); axh.set_xticks(np.arange(-0.5, len(ORD), 1), minor=True)
axh.tick_params(axis="x", which="major", length=2, width=0.6, pad=1)
axh.tick_params(axis="x", which="minor", length=0); axh.grid(False)
for sp in axh.spines.values(): sp.set_visible(False)
# 基因名（独立轴，完整显示）
axl.set_xlim(0, 1); axl.set_ylim(n-0.5, -0.5); axl.axis("off")
axg.set_xlim(0, 1); axg.set_ylim(n-0.5, -0.5); axg.axis("off")
for i, g in enumerate(P.gene):
    axl.text(1.0, i, g, ha="right", va="center", fontsize=5.3, fontstyle="italic", color=COL[P.grp[i]])
# 群间分隔线
b = 0
for g in ORD:
    k = int((P.grp == g).sum())
    if b > 0:
        for ax in (axl, axg, axs, axh):
            ax.axhline(b-0.5, color="white", lw=1.8)
    axg.text(1.0, b + k/2 - 0.5, g, ha="right", va="center", fontsize=6.2, weight="bold", color=COL[g])
    b += k
cb = fig.colorbar(im, cax=axc); cb.set_label("z-score", fontsize=6); cb.ax.tick_params(labelsize=5.5); cb.outline.set_visible(False)
fig.text(0.55, 0.975, "Marker genes", fontsize=8.5, ha="center")
bp = os.path.join(figdir, "Fig2_marker_heatmap")
for e, kw in (("_preview.png", dict(dpi=400)), (".pdf", {}), (".svg", {}), (".tiff", dict(dpi=1200, pil_kwargs={"compression":"tiff_lzw"}))):
    fig.savefig(bp + e, bbox_inches="tight", pad_inches=0.04, **kw)
plt.close()
for e in ("_preview.png",".pdf",".svg",".tiff"): print("  %-30s %6.2f MB" % (os.path.basename(bp+e), os.path.getsize(bp+e)/1e6))


