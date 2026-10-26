import sys, io, os
import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace"); pd.set_option("display.width", 250)
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
base = root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result"
A = sc.read_h5ad(out + r"\adata_raw.h5ad")
cl = pd.read_csv(base + r"\2.Cluster\2.cluster\Cells.cluster.list.xls", sep="\t"); cl["Cells"]=cl["Cells"].astype(str)
A = A[A.obs_names.isin(set(cl["Cells"]))].copy()
mp = cl.set_index("Cells")
A.obs["clu"] = mp.loc[A.obs_names,"Cluster"].astype(str).values
A.obs["group"] = mp.loc[A.obs_names,"Groups"].astype(str).values
A.obs["sample"] = mp.loc[A.obs_names,"Samples"].astype(str).values
# 给血细胞加上细分标签
Hh = sc.read_h5ad(out + r"\adata_hemocyte_harmony.h5ad")
A.obs["sub"] = "other"
A.obs.loc[Hh.obs_names, "sub"] = ("hemo_" + Hh.obs["sub_0.3"].astype(str)).values
A.layers["counts"] = A.layers["counts"].copy(); A.X = A.layers["counts"].copy()
sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
sc.pp.highly_variable_genes(A, n_top_genes=2000, flavor="seurat_v3", layer="counts")
G = A[:, A.var.highly_variable].copy(); G.X = G.layers["counts"].copy()
sc.pp.normalize_total(G, target_sum=1e4); sc.pp.log1p(G); sc.pp.scale(G, max_value=10)
sc.tl.pca(G, n_comps=50, svd_solver="arpack")
sc.external.pp.harmony_integrate(G, "sample", basis="X_pca", adjusted_basis="X_pca_harmony")
sc.pp.neighbors(G, use_rep="X_pca_harmony", n_neighbors=15, n_pcs=30)
# 全局 PAGA：非血细胞用原分群，血细胞用细分
G.obs["grp"] = np.where(G.obs["sub"] == "other", "C" + G.obs["clu"].astype(str), G.obs["sub"].astype(str))
G.obs["grp"] = G.obs["grp"].astype("category")
sc.tl.paga(G, groups="grp")
ct = pd.DataFrame(G.uns["paga"]["connectivities"].todense(), index=G.obs["grp"].cat.categories, columns=G.obs["grp"].cat.categories)
print("=== 全局 PAGA 连通性（组名 × 组名）==="); print(ct.round(3).to_string())
ct.to_csv(out + r"\paga_global_connectivity.csv")
if "hemo_2" in ct.index:
    print("\n=== 增殖祖细胞(hemo_2) 与其他组的连接强度（排序）===")
    print(ct.loc["hemo_2"].drop("hemo_2").sort_values(ascending=False).round(3).to_string())
print("\n=== hemo_2 与各组的最短路径（PAGA 图，越小越近）===")
import networkx as nx
Gg = nx.from_pandas_adjacency(ct); Gg.remove_edges_from([(u,u) for u in Gg.nodes])
import numpy as np
for t in ct.index:
    try: print("   hemo_2 -> %-10s 路径长度 %.0f" % (t, nx.shortest_path_length(Gg, "hemo_2", t)))
    except Exception: print("   hemo_2 -> %-10s 不连通" % t)
sc.tl.umap(G)
sc.settings.figdir = out + r"\figures_traj"
sc.pl.paga(G, color="grp", show=False, save="_paga_global.png", node_size_scale=2.0, edge_width_scale=3.0)
sc.pl.umap(G, color=["grp","sample"], show=False, save="_global.png", ncols=2)
G.write(out + r"\adata_global_integrated.h5ad")
print("\nfigures:", sorted(os.listdir(out + r"\figures_traj")))
