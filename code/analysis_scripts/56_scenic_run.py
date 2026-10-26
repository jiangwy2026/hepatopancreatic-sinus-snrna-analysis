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
from dask.distributed import Client, LocalCluster
cluster = LocalCluster(n_workers=8, threads_per_worker=1, processes=False, diagnostics_port=None)
client = Client(cluster)
print("dask 线程模式已启动", flush=True)
from arboreto.algo import grnboost2
t0 = time.time()
adj = grnboost2(expr, tf_names=use_tfs, client_or_address=client, seed=0, verbose=False)
adj.to_csv(out + r"\scenic_adjacencies.csv", index=False)
print("GRN: %d 条关系 (%.0fs)" % (len(adj), time.time()-t0), flush=True)
from pyscenic.utils import modules_from_adjacencies
mods = list(modules_from_adjacencies(adj, expr, rho_mask_dropouts=False))
print("模块: %d" % len(mods), flush=True)
from ctxcore.rnkdb import FeatherRankingDatabase as RDB
from pyscenic.prune import prune2df, df2regulons
rdb = RDB(os.path.join(DB, "dm6_v10_clust.genes_vs_motifs.rankings.feather"), "dm6")
motif_ann = pd.read_csv(os.path.join(DB, "motifs-v9-nr.flybase.tbl"), sep="\t", low_memory=False)
t1 = time.time()
try:
    df = prune2df([rdb], mods, motif_ann, num_workers=1)
except Exception as e:
    print("prune 失败:", str(e)[:150]); raise
regs = df2regulons(df)
print("regulon: %d 个 (%.0fs)" % (len(regs), time.time()-t1), flush=True)
import pickle; pickle.dump(regs, open(out + r"\scenic_regulons.pkl","wb"))
from pyscenic.aucell import aucell
auc = aucell(expr, regs, num_workers=1, seed=0)
auc.to_csv(out + r"\scenic_regulon_auc.csv")
print("AUCell:", auc.shape, flush=True)
df.to_csv(out + r"\scenic_motif_enrichment.csv")
client.close(); cluster.close()
