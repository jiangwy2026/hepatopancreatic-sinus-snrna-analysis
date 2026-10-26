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
P = pd.read_csv(out + r"\marker_heatmap_selected.csv")
symv = A.var["sym"].astype(str).values; gidv = A.var_names.astype(str)
sym2i = {}; sym2g = {}
for i, (s, g) in enumerate(zip(symv, gidv)): sym2i.setdefault(s, i); sym2g.setdefault(s, g)
P["gid"] = P.gene.map(sym2g)
# 每群抽样细胞（大群最多 400 个，保持图幅可读）
rng = np.random.default_rng(0); cell_idx, bounds, marks = [], [], []
for g in ORD:
    idx = np.where((A.obs["ct"] == g).values)[0]
    if len(idx) > 200: idx = rng.choice(idx, 200, replace=False)
    idx = np.sort(idx); bounds.append(len(cell_idx)); cell_idx += list(idx)
    marks.append(len(cell_idx))
cell_idx = np.array(cell_idx)
E = A[cell_idx][:, [sym2i[g] for g in P.gene]].X.toarray().astype(np.float32)   # cells x genes (log-normalized)
# 参考图风格：按基因自身表达量归一化 -> 不表达=黑、高表达=黄（不做跨细胞 z-score）
den = np.percentile(E, 99, axis=0) + 1e-9
Z = np.clip(E / den, 0, 1).T                                                     # genes x cells, 0~1
print("热图矩阵: %d 基因 x %d 细胞 (每群 ≤200)" % Z.shape)
fig = plt.figure(figsize=(8.2, 5.4))
gs = fig.add_gridspec(2, 3, height_ratios=[0.055, 1.0], width_ratios=[0.26, 1.0, 0.030],
                      hspace=0.012, wspace=0.03, left=0.02, right=0.97, top=0.965, bottom=0.035)
axb = fig.add_subplot(gs[0, 1]); axh = fig.add_subplot(gs[1, 1])
axl = fig.add_subplot(gs[1, 0]); axc = fig.add_subplot(gs[1, 2])
nG, nC = Z.shape
# 顶部群色条
bars = np.zeros((1, nC))
for i, g in enumerate(ORD):
    a = 0 if i == 0 else marks[i-1]; bars[0, a:marks[i]] = i
axb.imshow(bars, aspect="auto", cmap=ListedColormap([COL[g] for g in ORD]),
           norm=BoundaryNorm(np.arange(-0.5, len(ORD)+0.5), len(ORD)), extent=(0, nC, 0, 1))
axb.set_xticks([]); axb.set_yticks([])
for i, g in enumerate(ORD):
    a = 0 if i == 0 else marks[i-1]; b = marks[i]
    axb.text((a+b)/2, 0.5, g, ha="center", va="center", fontsize=6, weight="bold", color="white", clip_on=False)
    if i > 0: axb.axvline(a, color="white", lw=1.2)
for sp in axb.spines.values(): sp.set_visible(False)
# 热图
im = axh.imshow(Z, aspect="auto", cmap="magma", vmin=0, vmax=1, interpolation="nearest")
axh.set_xticks([]); axh.set_yticks([])
for sp in axh.spines.values(): sp.set_visible(False)
for i in range(1, len(ORD)):
    a = marks[i-1]; axh.axvline(a, color="white", lw=1.0)
b = 0
for g in ORD:
    k = int((P.grp == g).sum())
    if b > 0: axh.axhline(b-0.5, color="white", lw=1.0)
    b += k
# 左侧基因名
axl.set_xlim(0, 1); axl.set_ylim(nG-0.5, -0.5); axl.axis("off")
for i, r in enumerate(P.itertuples()):
    axl.text(1.0, i, "%s (%s)" % (r.gene, r.gid) if pd.notna(r.gid) else r.gene,
             ha="right", va="center", fontsize=4.6, color="black")
b = 0
for g in ORD:
    k = int((P.grp == g).sum())
    if b > 0: axl.axhline(b-0.5, color="white", lw=1.0)
    b += k
cb = fig.colorbar(im, cax=axc); cb.set_label("Relative expression", fontsize=6); cb.ax.tick_params(labelsize=5.5); cb.outline.set_visible(False)
bp = os.path.join(figdir, "Fig2_marker_heatmap_singlecell")
for e, kw in (("_preview.png", dict(dpi=300)), (".pdf", {}), (".svg", {}), (".tiff", dict(dpi=1200, pil_kwargs={"compression":"tiff_lzw"}))):
    fig.savefig(bp + e, bbox_inches="tight", pad_inches=0.03, **kw)
plt.close()
for e in ("_preview.png",".pdf",".svg",".tiff"): print("  %-38s %6.2f MB" % (os.path.basename(bp+e), os.path.getsize(bp+e)/1e6))




