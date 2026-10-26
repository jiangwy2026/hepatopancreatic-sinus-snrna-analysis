import sys, io, os
import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"; fig = out + r"\figures"; os.makedirs(fig, exist_ok=True)
B = sc.read_h5ad(out + r"\adata_hvg_processed.h5ad")
sc.settings.figdir = fig
sc.pl.umap(B, color="vendor_cluster", legend_loc="on data", frameon=True, title="Vendor clusters (0-9)", save="_vendor_cluster.png", show=False, palette="tab10")
sc.pl.umap(B, color="group", frameon=True, title="Group", save="_group.png", show=False)
sc.pl.umap(B, color="sample", frameon=True, title="Sample", save="_sample.png", show=False)
print("figures:", sorted(os.listdir(fig)))

# marker panel: 用全基因表达谱（log-normalized）
A = sc.read_h5ad(out + r"\adata_raw.h5ad")
cl = pd.read_csv(root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\2.Cluster\2.cluster\Cells.cluster.list.xls", sep="\t")
cl["Cells"] = cl["Cells"].astype(str)
A = A[A.obs_names.isin(set(cl["Cells"]))].copy()
A.obs["clu"] = cl.set_index("Cells").loc[A.obs_names, "Cluster"].astype(str).values
A.X = A.layers["counts"].copy(); sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
panel = {
 "hemocyte": ["A2M","A2ML1","Pxt","PPAF1","PPAF2","LOC113823945","LOC113801825","LOC113816267","PEN-3","LOC113800363","LOC113802172","AGAP007347","LOC113820206","LOC113813598","LOC113815775"],
 "proliferation": ["ccnd2","Pka-C1","Patronin"],
 "hepatopancreas": ["Prss1","Prss2","Try3","Jon99Cii","Cpa2","Chia","LCP2","LIPF","ACE","Mal-A2","Oacyl","hdc"],
 "cuticle": ["LOC113809538","LOC113809539","LOC113809540","LOC113819135","LOC113810024","LOC113810025"],
 "neural": ["RASGRF1","ADCY5","LOC113803948","KCNQ5","Nrg","Fas2","Grik2","slo"],
}
sym2id = {}
for i, s in zip(A.var_names, A.var["symbol"]):
    sym2id.setdefault(s, i)
ok, miss = {}, []
for grp, genes in panel.items():
    for g in genes:
        if g in sym2id and g not in ok: ok[g] = sym2id[g]
        elif g not in sym2id: miss.append(g)
print("\n未找到的基因:", miss)
ids = list(ok.values()); names = list(ok.keys())
sub = A[:, A.var_names.isin(ids)].copy()
df = pd.DataFrame(np.asarray(sub.X.todense()), columns=list(sub.var["symbol"].values))
df["clu"] = A.obs["clu"].values
M = df.groupby("clu").mean().T
M = M.loc[[n for n in names if n in M.index]]
M.to_csv(out + r"\marker_panel_cluster_means.csv")
pd.set_option("display.width", 250)
print("\n=== marker 面板：各群平均表达（log-normalized）===")
print(M.round(2).to_string())
