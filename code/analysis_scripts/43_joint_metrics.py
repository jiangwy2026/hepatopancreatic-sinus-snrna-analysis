import sys, io, os
import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, scanpy as sc
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 260); sc.settings.verbosity = 0
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
base = root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result"
H = sc.read_h5ad(out + r"\joint_mine_elife.h5ad")
E = sc.read_h5ad(out + r"\elife_adata.h5ad")
cl = pd.read_csv(base + r"\2.Cluster\2.cluster\Cells.cluster.list.xls", sep="\t"); cl["Cells"]=cl["Cells"].astype(str)
myc_series = cl.set_index("Cells")["Cluster"].astype(str)
H.obs["elife_cluster"] = E.obs["elife_cluster"].astype(str).reindex(H.obs_names).values
H.obs["my_cluster"] = myc_series.reindex(H.obs_names).values
H.obs["my_cluster"] = H.obs["my_cluster"].fillna("NA").astype(str)
mine = (H.obs.dataset == "Mine").values; el = ~mine
elc = H.obs["elife_cluster"].values; myc = H.obs["my_cluster"].values
W = H.obsp["connectivities"].tocsr()
print("我的细胞 %d | eLife 细胞 %d" % (mine.sum(), el.sum()))
rows = []
for i in np.where(mine)[0]:
    nb = W.indices[W.indptr[i]:W.indptr[i+1]]; enb = nb[el[nb]]
    vc = pd.Series(elc[enb]).value_counts() if len(enb) else None
    rows.append({"my_cluster": myc[i], "n_nb": len(nb), "n_elife": len(enb),
                 "frac_elife": len(enb)/max(len(nb),1),
                 "top_elife": (vc.index[0] if vc is not None and len(vc) else "NA"),
                 "top_share": (vc.iloc[0]/vc.sum() if vc is not None and len(vc) else np.nan)})
R = pd.DataFrame(rows)
g = R.groupby("my_cluster").agg(n=("frac_elife","size"), 邻居含eLife比例=("frac_elife","mean"),
        最常邻接eLife群=("top_elife", lambda s: s.mode().iloc[0] if len(s.mode()) else "NA"),
        该群占邻居比=("top_share","mean"))
print("\n=== 我的每个群：邻居里 eLife 细胞比例 / 最常邻接的 eLife 群 ===")
print(g.round(3).to_string())
rows2 = []
for i in np.where(el)[0]:
    nb = W.indices[W.indptr[i]:W.indptr[i+1]]; mnb = nb[mine[nb]]
    vc = pd.Series(myc[mnb]).value_counts() if len(mnb) else None
    rows2.append({"elife_cluster": elc[i], "frac_mine": len(mnb)/max(len(nb),1),
                  "top_mine": (vc.index[0] if vc is not None and len(vc) else "NA"),
                  "top_share": (vc.iloc[0]/vc.sum() if vc is not None and len(vc) else np.nan)})
R2 = pd.DataFrame(rows2)
g2 = R2.groupby("elife_cluster").agg(n=("frac_mine","size"), 邻居含我的比例=("frac_mine","mean"),
        最常邻接我的群=("top_mine", lambda s: s.mode().iloc[0] if len(s.mode()) else "NA"),
        该群占邻居比=("top_share","mean"))
print("\n=== 反过来：eLife 各群邻接我的哪个群最多 ==="); print(g2.round(3).to_string())
R.to_csv(out + r"\joint_mine_neighbors.csv", index=False); R2.to_csv(out + r"\joint_elife_neighbors.csv", index=False)
um = H.obsm["X_umap"]
print("\n=== 联合 UMAP 质心位置 ===")
print("  我的各群:"); print(pd.DataFrame({"c": myc[mine], "x": um[mine,0], "y": um[mine,1]}).groupby("c").mean().round(1).to_string())
print("  eLife 各群:"); print(pd.DataFrame({"c": elc[el], "x": um[el,0], "y": um[el,1]}).groupby("c").mean().round(1).to_string())
