import sys, io, os, time
import numpy as np, pandas as pd, scanpy as sc
# numpy 2.x 移除了 np.object 等别名，pyscenic 0.12.1 仍在用
for _a, _b in (("object", object), ("bool", bool), ("int", int), ("float", float), ("str", str)):
    if not hasattr(np, _a): setattr(np, _a, _b)
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
pd.set_option("display.width", 250); sc.settings.verbosity = 0
root = r"D:\桌面\single cell"; out = root + r"\analysis_output"; DB = out + r"\scenic_db"
# --- 打补丁：pyscenic 0.12.1 在 pandas 3.0 下 groupby().apply() 会丢掉分组列 ---
import pyscenic.utils as u
from ctxcore.genesig import Regulon
def modules4top_factors_fixed(adjacencies, n, context=frozenset()):
    df = adjacencies.sort_values(u.COLUMN_NAME_WEIGHT, ascending=False)
    df = df.groupby(u.COLUMN_NAME_TARGET, sort=False, as_index=False).head(n)
    for tf_name, df_grp in df.groupby(by=u.COLUMN_NAME_TF):
        if len(df_grp) > 0:
            yield Regulon(name=tf_name, context=frozenset(["top%dperTarget" % n]).union(context),
                          transcription_factor=tf_name,
                          gene2weight=list(zip(df_grp[u.COLUMN_NAME_TARGET].values,
                                               df_grp[u.COLUMN_NAME_WEIGHT].values)),
                          gene2occurrence=[])
u.modules4top_factors = modules4top_factors_fixed
print("补丁已应用", flush=True)
H = sc.read_h5ad(out + r"\adata_hemocyte_traj.h5ad")
H.X = H.layers["counts"].copy(); sc.pp.normalize_total(H, target_sum=1e4); sc.pp.log1p(H)
dbgenes = set(pd.read_csv(out + r"\scenic_gene_overlap.csv", header=None)[0].astype(str))
tfset = set(pd.read_csv(out + r"\scenic_TF_list.csv", header=None)[0].astype(str))
sym2i = {}
for i, s in zip(H.var_names, H.var["sym"].astype(str)):
    if s in dbgenes: sym2i.setdefault(s, i)
genes = sorted(sym2i)
expr = pd.DataFrame(H[:, [sym2i[g] for g in genes]].X.toarray(), index=H.obs_names, columns=genes)
adj = pd.read_csv(out + r"\scenic_adjacencies.csv")
adj = adj[adj.TF.isin(expr.columns) & adj.target.isin(expr.columns)].copy()
print("adjacencies:", adj.shape, flush=True)
mods = list(u.modules_from_adjacencies(adj, expr, rho_mask_dropouts=False))
print("模块: %d" % len(mods), flush=True)
from ctxcore.rnkdb import FeatherRankingDatabase as RDB
from pyscenic.prune import prune2df, df2regulons
rdb = RDB(os.path.join(DB, "dm6_v10_clust.genes_vs_motifs.rankings.feather"), "dm6")
ann = pd.read_csv(os.path.join(DB, "motifs-v9-nr.flybase.tbl"), sep="\t", low_memory=False)
t0 = time.time(); df = prune2df([rdb], mods, ann, num_workers=1)
regs = df2regulons(df)
print("regulon: %d 个 (%.0fs)" % (len(regs), time.time()-t0), flush=True)
print("regulon 列表:", [r.name for r in regs][:40], flush=True)
df.to_csv(out + r"\scenic_motif_enrichment.csv")
from pyscenic.aucell import aucell
auc = aucell(expr, regs, num_workers=1, seed=0)
auc.to_csv(out + r"\scenic_regulon_auc.csv")
print("AUCell:", auc.shape, flush=True)
import pickle; pickle.dump([(r.name, r.transcription_factor, len(r)) for r in regs], open(out + r"\scenic_regulon_info.pkl","wb"))


