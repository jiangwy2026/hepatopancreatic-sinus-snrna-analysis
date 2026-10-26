import sys, io, os
import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, scanpy as sc
from scipy.stats import spearmanr
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace"); pd.set_option("display.width", 250)
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"; fig = out + r"\figures_traj"; os.makedirs(fig, exist_ok=True)
S_G = ['PCNA','MCM6','CDC6','ORC1','ORC4','ORC5','POLA2','PRIM1','FEN1','TYMS','DHFR','HIRA','CHAF1A']
G2M_G = ['CCNB3','AURKA','TOP2A','TPX2','KIF20A','KIF23','RACGAP1','NDC80','INCENP','PRC1','ASPM','Pka-C1','Patronin']
Hh = sc.read_h5ad(out + r"\adata_hemocyte_harmony.h5ad")
A = sc.read_h5ad(out + r"\adata_raw.h5ad"); A = A[list(Hh.obs_names)].copy()
A.obsm["X_pca_harmony"] = Hh.obsm["X_pca_harmony"]
A.obs["sub"] = Hh.obs["sub_0.3"].astype(str).values
A.obs["group"] = Hh.obs["group"].astype(str).values
A.X = A.layers["counts"].copy(); sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
sym2id = {}
for i, s in zip(A.var_names, A.var["symbol"]): sym2id.setdefault(str(s), i)
sc.tl.score_genes_cell_cycle(A, s_genes=[sym2id[g] for g in S_G], g2m_genes=[sym2id[g] for g in G2M_G])
sc.pp.neighbors(A, use_rep="X_pca_harmony", n_neighbors=15, n_pcs=30)
sc.tl.paga(A, groups="sub")
sc.tl.diffmap(A, n_comps=15)
sc_ = A.obs.groupby("sub")[["S_score","G2M_score"]].mean(); sc_["tot"] = sc_.sum(axis=1); rootlab = sc_["tot"].idxmax()
idx = np.where((A.obs["sub"] == rootlab).values)[0]
iroot = int(idx[np.argmax(A.obs["G2M_score"].values[idx])]); A.uns["iroot"] = iroot
sc.tl.dpt(A)
sc.tl.umap(A)
print("root 亚群 =", rootlab, "| root 细胞 =", A.obs_names[iroot])
print("\n=== 各亚群 pseudotime ==="); print(A.obs.groupby("sub")["dpt_pseudotime"].agg(["mean","median","min","max"]).round(3).to_string())
print("\n=== blank vs treated pseudotime ===")
for g in ["blank","treated"]:
    v = A.obs.loc[A.obs.group==g,"dpt_pseudotime"]
    print("  %-8s 中位=%.3f | 早期<0.2: %.1f%% | 0.2-0.5: %.1f%% | >0.5: %.1f%%" % (g, v.median(), (v<0.2).mean()*100, ((v>=0.2)&(v<=0.5)).mean()*100, (v>0.5).mean()*100))
print("\n=== 各亚群 x 组的 pseudotime 中位数 ===")
print(A.obs.pivot_table(index="sub", columns="group", values="dpt_pseudotime", aggfunc="median").round(3).to_string())
X = A.X if not hasattr(A.X,"todense") else np.asarray(A.X.todense())
hv = np.argsort(-X.var(0))[:3000]; pt = A.obs["dpt_pseudotime"].values
rng = np.random.default_rng(0); s = rng.choice(len(pt), size=min(4000,len(pt)), replace=False)
rho, pv = spearmanr(X[s][:, hv], pt[s])
res = pd.DataFrame({"gene": A.var["symbol"].values[hv].astype(str), "rho": rho, "p": pv}).dropna().sort_values("rho")
res.to_csv(out + r"\genes_vs_pseudotime.csv", index=False)
ann = {}
with open(root + r"\all.genes.marker.stat.xls.xls", encoding="utf-8", errors="replace") as fh:
    fh.readline()
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 24: ann[p[1]] = p[14][:55]
print("\n=== pseudotime 下降基因（祖细胞特征）===")
for _, r in res.head(18).iterrows(): print("   %-20s rho=%+.2f | %s" % (r.gene, r.rho, ann.get(r.gene,"")))
print("\n=== pseudotime 上升基因（成熟细胞特征）===")
for _, r in res.tail(18).iloc[::-1].iterrows(): print("   %-20s rho=%+.2f | %s" % (r.gene, r.rho, ann.get(r.gene,"")))
sc.settings.figdir = fig
sc.pl.paga(A, color="sub", show=False, save="_paga.png", title="PAGA (hemocytes)")
sc.pl.umap(A, color=["sub","phase","dpt_pseudotime","group"], show=False, save="_traj.png", ncols=2)
sc.pl.draw_graph(A, color=["sub","dpt_pseudotime"], show=False, save="_graph.png")
A.write(out + r"\adata_hemocyte_dpt.h5ad")
print("\nfigures:", sorted(os.listdir(fig)))
