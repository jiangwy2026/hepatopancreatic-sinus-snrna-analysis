import sys, io, os, time
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"
EL = out + r"\elife_data\GDR21050291_std_1-Penaeus_vannamei-scRNAseq-result"
myfeat = pd.read_csv(root + r"\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\ADRL26080060_std_1-Litopenaeus_vannamei-scRNAseq-result\1.Expression\expressions\blank\features.tsv.gz", sep="\t", header=None, names=["gid","sym","type"])
eg = pd.read_csv(EL + r"\1.Expression\expressions\Ctrl\genes.tsv", sep="\t", header=None)
print("=== 基因 ID 对比 ===")
print("  eLife 基因数: %d | 你的基因数: %d" % (len(eg), len(myfeat)))
print("  eLife 前3:", eg[0].head(3).tolist())
print("  你的 前3:", myfeat["gid"].head(3).tolist())
inter = set(eg[0].astype(str)) & set(myfeat["gid"].astype(str))
print("  完全相同的 ID 交集: %d (%.1f%%)" % (len(inter), len(inter)/len(myfeat)*100))
print("\n=== 载入 eLife 3 个样本 ===")
parts = []
for s in ("Ctrl","rEGFP","rCREG"):
    t0 = time.time(); base = EL + r"\1.Expression\expressions" + "\\" + s
    bars = [l.strip() for l in open(base + r"\barcodes.tsv", encoding="utf-8") if l.strip()]
    with open(base + r"\matrix.mtx") as f:
        for line in f:
            if line.startswith("%"): continue
            dims = line.split(); break
    d = pd.read_csv(base + r"\matrix.mtx", sep=" ", skiprows=3, header=None, names=["g","c","x"],
                    dtype={"g":np.int32,"c":np.int32,"x":np.int32}, engine="c")
    M = sp.csr_matrix((d["x"].to_numpy(np.float32), (d["c"].to_numpy()-1, d["g"].to_numpy()-1)),
                      shape=(int(dims[1]), int(dims[0])))
    A = ad.AnnData(X=M)
    A.var_names = eg[0].astype(str).values
    A.obs_names = [s + "_" + b.rsplit("-", 1)[0] for b in bars]
    A.obs["sample"] = s
    print("  %-6s %5d 细胞 x %5d 基因 | nnz=%d | 中位基因数 %.0f | %.0fs" %
          (s, A.n_obs, A.n_vars, M.nnz, np.median(np.asarray((M>0).sum(1)).ravel()), time.time()-t0), flush=True)
    parts.append(A); del d, M
E = ad.concat(parts, join="outer", index_unique=None)
del parts
E.var["sym"] = myfeat.set_index("gid").reindex(E.var_names)["sym"].values
cl = pd.read_csv(EL + r"\2.Cluster\2.cluster\Cells.cluster.list.xls", sep="\t")
cl["Cells"] = cl["Cells"].astype(str)
E.obs["elife_cluster"] = cl.set_index("Cells").reindex(E.obs_names)["Cluster"].astype("Int64").astype(str).values
cc = pd.read_csv(EL + r"\5.CellCycle\CellCycle.annot.xls", sep="\t")
cc["Cells"] = cc["Cells"].astype(str)
E.obs["phase"] = cc.set_index("Cells").reindex(E.obs_names)["Phase"].values
E.layers["counts"] = E.X.copy()
print("\neLife 合并:", E.shape, "| 分群数:", E.obs.elife_cluster.nunique())
print("分群分布:"); print(E.obs.elife_cluster.value_counts().sort_index(key=lambda x: x.astype(int)).to_string())
print("相位分布:"); print(E.obs.phase.value_counts().to_string())
E.write(out + r"\elife_adata.h5ad")
print("\nsaved elife_adata.h5ad")
