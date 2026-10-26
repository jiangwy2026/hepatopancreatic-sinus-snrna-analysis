import sys, io, os, traceback
import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, scanpy as sc
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace"); pd.set_option("display.width", 250)
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"; fig = out + r"\figures_traj"; os.makedirs(fig, exist_ok=True)
A = sc.read_h5ad(out + r"\adata_raw.h5ad")
D = pd.read_csv(out + r"\fate_probabilities.csv") if os.path.exists(out + r"\fate_probabilities.csv") else None
import subprocess
sub = pd.read_csv(out + r"\subcluster_markers_final.csv") if False else None
# 重新构造血细胞对象（含 harmony 坐标与亚群标签）
Hh = sc.read_h5ad(out + r"\adata_hemocyte_harmony.h5ad")
A = A[list(Hh.obs_names)].copy()
A.obsm["X_pca_harmony"] = Hh.obsm["X_pca_harmony"]
A.obs["sub"] = Hh.obs["sub_0.3"].astype(str).values
A.obs["group"] = Hh.obs["group"].astype(str).values
A.X = A.layers["counts"].copy(); sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
S_G = ['PCNA','MCM6','CDC6','ORC1','ORC4','ORC5','POLA2','PRIM1','FEN1','TYMS','DHFR','HIRA','CHAF1A']
G2M_G = ['CCNB3','AURKA','TOP2A','TPX2','KIF20A','KIF23','RACGAP1','NDC80','INCENP','PRC1','ASPM','Pka-C1','Patronin']
sym2id = {}
for i, s in zip(A.var_names, A.var["symbol"]): sym2id.setdefault(str(s), i)
sc.tl.score_genes_cell_cycle(A, s_genes=[sym2id[g] for g in S_G], g2m_genes=[sym2id[g] for g in G2M_G])
Al = A[A.obs["sub"].isin(["0","1","2"])].copy()
sc.pp.neighbors(Al, use_rep="X_pca_harmony", n_neighbors=15, n_pcs=30)
sc.tl.paga(Al, groups="sub")
sc.tl.diffmap(Al, n_comps=15)
sc_ = Al.obs.groupby("sub")[["S_score","G2M_score"]].mean(); sc_["tot"]=sc_.sum(axis=1); rootlab = sc_["tot"].idxmax()
idx = np.where((Al.obs["sub"]==rootlab).values)[0]; iroot = int(idx[np.argmax(Al.obs["G2M_score"].values[idx])])
Al.uns["iroot"] = iroot; sc.tl.dpt(Al); sc.tl.umap(Al)
print("root 亚群 =", rootlab, "| pseudotime 各亚群中位:"); print(Al.obs.groupby("sub")["dpt_pseudotime"].median().round(3).to_string())

# ---- 基因随 pseudotime 变化（秩相关，向量化）----
from scipy.stats import rankdata
X = Al.X if not hasattr(Al.X,"todense") else np.asarray(Al.X.todense())
pt = Al.obs["dpt_pseudotime"].values
keep = pt > 0
Xo = X[keep]; pto = pt[keep]
rng = np.random.default_rng(0); s_idx = rng.choice(Xo.shape[0], size=min(5000, Xo.shape[0]), replace=False)
Xs = Xo[s_idx]; pts = pto[s_idx]
ry = rankdata(pts); ry = (ry - ry.mean()) / (ry.std() + 1e-12)
var = Xs.var(0); hv = np.argsort(-var)[:4000]
rx = rankdata(Xs[:, hv], axis=0); rx = (rx - rx.mean(0)) / (rx.std(0) + 1e-12)
rho = (rx * ry[:, None]).mean(0)
res = pd.DataFrame({"gene": Al.var["symbol"].values[hv].astype(str), "rho": rho}).sort_values("rho")
res.to_csv(out + r"\genes_vs_pseudotime.csv", index=False)
ann = {}
with open(root + r"\all.genes.marker.stat.xls.xls", encoding="utf-8", errors="replace") as fh:
    fh.readline()
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 24: ann[p[1]] = p[14][:52]
print("\n=== 随拟时序下降（祖细胞/增殖特征）top15 ===")
for _, r in res.head(15).iterrows(): print("   %-20s rho=%+.2f | %s" % (r.gene, r.rho, ann.get(r.gene,"")))
print("\n=== 随拟时序上升（成熟/分化特征）top15 ===")
for _, r in res.tail(15).iloc[::-1].iterrows(): print("   %-20s rho=%+.2f | %s" % (r.gene, r.rho, ann.get(r.gene,"")))
Al.write(out + r"\adata_hemocyte_lineage.h5ad")
sc.settings.figdir = fig
sc.pl.paga(Al, color="sub", show=False, save="_paga.png")
sc.pl.umap(Al, color=["sub","phase","dpt_pseudotime"], show=False, save="_traj.png", ncols=3)
sc.pl.umap(A, color=["sub","phase","group"], show=False, save="_all6.png", ncols=3)
print("\nfigures:", sorted(os.listdir(fig)))
