import sys, io, os, time
import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, scanpy as sc
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 250); sc.settings.verbosity = 0
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"; DB = out + r"\scenic_db"
H = sc.read_h5ad(out + r"\adata_hemocyte_traj.h5ad")
H.X = H.layers["counts"].copy(); sc.pp.normalize_total(H, target_sum=1e4); sc.pp.log1p(H)
dbgenes = set(pd.read_csv(out + r"\scenic_gene_overlap.csv", header=None)[0].astype(str))
tfs = set(pd.read_csv(out + r"\scenic_TF_list.csv", header=None)[0].astype(str))
sym2i = {}
for i, s in zip(H.var_names, H.var["sym"].astype(str)):
    if s in dbgenes: sym2i.setdefault(s, i)
genes = sorted(sym2i); use_tfs = sorted(set(genes) & tfs)
expr = pd.DataFrame(H[:, [sym2i[g] for g in genes]].X.toarray(), index=H.obs_names, columns=genes)
print("输入: %d 细胞 x %d 基因 | TF %d" % (expr.shape[0], expr.shape[1], len(use_tfs)), flush=True)
adj = None
# 方案 A: 自己建 dask 客户端（不带 diagnostics_port）
try:
    from dask.distributed import Client, LocalCluster
    cl = LocalCluster(n_workers=8, threads_per_worker=1, processes=False)
    client = Client(cl)
    from arboreto.algo import grnboost2
    t0 = time.time()
    adj = grnboost2(expr, tf_names=use_tfs, client_or_address=client, seed=0, verbose=False)
    print("方案A(GRNBoost2) 成功: %d 条关系 (%.0fs)" % (len(adj), time.time()-t0), flush=True)
    client.close(); cl.close()
except Exception as e:
    print("方案A 失败:", str(e)[:140], flush=True)
# 方案 B: 手动 GBM（等价于 GRNBoost2 的核心算法）
if adj is None:
    from sklearn.ensemble import GradientBoostingRegressor
    t0 = time.time(); rows = []
    Xall = expr.values.astype(np.float32); names = list(expr.columns)
    idx_all = np.arange(len(names))
    for k, tf in enumerate(use_tfs):
        j = names.index(tf)
        y = Xall[:, j]; X = np.delete(Xall, j, axis=1); tn = [n for n in names if n != tf]
        m = GradientBoostingRegressor(n_estimators=100, learning_rate=0.01, max_depth=3, subsample=0.9, random_state=0)
        m.fit(X, y)
        imp = m.feature_importances_
        top = np.argsort(-imp)[:50]
        for t in top:
            if imp[t] > 0: rows.append((tf, tn[t], float(imp[t])))
        if (k+1) % 50 == 0: print("  已处理 %d/%d TF (%.0fs)" % (k+1, len(use_tfs), time.time()-t0), flush=True)
    adj = pd.DataFrame(rows, columns=["TF","target","importance"])
    print("方案B(GBM) 完成: %d 条关系 (%.0fs)" % (len(adj), time.time()-t0), flush=True)
adj.to_csv(out + r"\scenic_adjacencies.csv", index=False)
from pyscenic.utils import modules_from_adjacencies
mods = list(modules_from_adjacencies(adj, expr, rho_mask_dropouts=False))
print("模块: %d" % len(mods), flush=True)
from ctxcore.rnkdb import FeatherRankingDatabase as RDB
from pyscenic.prune import prune2df, df2regulons
rdb = RDB(os.path.join(DB, "dm6_v10_clust.genes_vs_motifs.rankings.feather"), "dm6")
motif_ann = pd.read_csv(os.path.join(DB, "motifs-v9-nr.flybase.tbl"), sep="\t", low_memory=False)
t1 = time.time()
df = prune2df([rdb], mods, motif_ann, num_workers=1)
regs = df2regulons(df)
print("regulon: %d 个 (%.0fs)" % (len(regs), time.time()-t1), flush=True)
import pickle; pickle.dump([r.name for r in regs], open(out + r"\scenic_regulon_names.pkl","wb"))
from pyscenic.aucell import aucell
auc = aucell(expr, regs, num_workers=1, seed=0)
auc.to_csv(out + r"\scenic_regulon_auc.csv")
df.to_csv(out + r"\scenic_motif_enrichment.csv")
print("AUCell:", auc.shape, flush=True)
print("regulon 列表(前30):", [r.name for r in regs][:30], flush=True)
