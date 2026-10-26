import sys, io
import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, scanpy as sc
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace"); pd.set_option("display.width", 250)
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
S_G = ['PCNA','MCM6','CDC6','ORC1','ORC4','ORC5','POLA2','PRIM1','FEN1','TYMS','DHFR','HIRA','CHAF1A']
G2M_G = ['CCNB3','AURKA','TOP2A','TPX2','KIF20A','KIF23','RACGAP1','NDC80','INCENP','PRC1','ASPM','Pka-C1','Patronin']

Hh = sc.read_h5ad(out + r"\adata_hemocyte_harmony.h5ad")
A  = sc.read_h5ad(out + r"\adata_raw.h5ad")
A = A[list(Hh.obs_names)].copy()
A.obsm["X_pca_harmony"] = Hh.obsm["X_pca_harmony"]
A.obs["sub"]  = Hh.obs["sub_0.3"].astype(str).values
A.obs["group"]= Hh.obs["group"].astype(str).values
A.X = A.layers["counts"].copy(); sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
sym2id = {}
for i, s in zip(A.var_names, A.var["symbol"]): sym2id.setdefault(str(s), i)
sc.tl.score_genes_cell_cycle(A, s_genes=[sym2id[g] for g in S_G], g2m_genes=[sym2id[g] for g in G2M_G])
sc.pp.neighbors(A, use_rep="X_pca_harmony", n_neighbors=15, n_pcs=30, key_added="hn")
lab = "sub"
print("=== 各亚群细胞周期分布（%）===")
ph = pd.crosstab(A.obs[lab], A.obs["phase"], normalize="index").mul(100).round(1); ph["n"] = A.obs[lab].value_counts()
print(ph.to_string())
print("\n=== 各亚群 S/G2M 打分均值 ===")
print(A.obs.groupby(lab)[["S_score","G2M_score"]].mean().round(3).to_string())
sc.tl.paga(A, groups=lab, neighbors_key="hn")
ct = pd.DataFrame(A.uns["paga"]["connectivities"].todense(), index=A.obs[lab].cat.categories, columns=A.obs[lab].cat.categories)
print("\n=== PAGA 连通性（行=起点 → 列=终点）==="); print(ct.round(3).to_string()); ct.to_csv(out + r"\paga_connectivity_hemocyte.csv")
score = A.obs.groupby(lab)[["S_score","G2M_score"]].mean().sum(1).sort_values(ascending=False)
rootlab = score.index[0]; print("\n周期性最强亚群 =", rootlab, dict(score.round(3)))
sc.tl.diffmap(A, n_comps=15, neighbors_key="hn")
idx = np.where((A.obs[lab] == rootlab).values)[0]
iroot = int(idx[np.argmax(A.obs["G2M_score"].values[idx])]); A.uns["iroot"] = iroot
sc.tl.dpt(A)
print("root:", A.obs_names[iroot], "亚群", A.obs[lab].iloc[iroot])
print("\n=== 各亚群 dpt_pseudotime ==="); print(A.obs.groupby(lab)["dpt_pseudotime"].agg(["mean","median"]).round(3).to_string())
print("\n=== blank vs treated 的 pseudotime ===")
for g in ["blank","treated"]:
    v = A.obs.loc[A.obs.group==g,"dpt_pseudotime"]
    print("  %-8s 25%%=%.3f 中位=%.3f 75%%=%.3f | 晚期(>0.5)占比=%.1f%% | 早期(<0.2)=%.1f%%" % (g, v.quantile(.25), v.median(), v.quantile(.75), (v>0.5).mean()*100, (v<0.2).mean()*100))
A.write(out + r"\adata_hemocyte_dpt.h5ad")
print("\nsaved adata_hemocyte_dpt.h5ad")
