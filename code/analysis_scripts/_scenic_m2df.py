# -*- coding: utf-8 -*-
import os, inspect, numpy as np, pandas as pd
np.object = object
import dask; dask.config.set({"dataframe.query-planning": False})
from pyscenic import prune as P
print("funcs:", [n for n in dir(P) if "df" in n or "module" in n.lower()])
print("modules2df:", inspect.signature(P.modules2df))
print("prune2df:", inspect.signature(P.prune2df))
OUT = r"D:\桌面\single cell\analysis_output"; DB = os.path.join(OUT, "scenic_db")
expr = pd.read_csv(r"E:\workspace\scenic_work\expr.csv.gz", index_col=0)
adj = pd.read_csv(os.path.join(OUT, "scenic_adjacencies.csv")); adj.columns = ["TF","target","importance"]
adj = adj[adj.TF.isin(expr.columns) & adj.target.isin(expr.columns)]
from pyscenic.utils import modules_from_adjacencies
mods = list(modules_from_adjacencies(adj, expr, rho_mask_dropouts=False))
print("modules:", len(mods), "| sizes:", [len(m) for m in mods[:10]])
from ctxcore.rnkdb import FeatherRankingDatabase as RDB
db = RDB(os.path.join(DB, "dm6_v10_clust.genes_vs_motifs.rankings.feather"), "dm6")
ann = pd.read_csv(os.path.join(DB, "motifs-v10nr_clust-nr.flybase-m0.001-o0.0.tbl"), sep="\t", low_memory=False)
sub = mods[:5]
for nes in (3.0, 0.0):
    try:
        df = P.modules2df(db, sub, ann, rank_threshold=1500, auc_threshold=0.05, nes_threshold=nes)
    except TypeError as e:
        print("signature issue:", repr(e)[:200])
        df = P.modules2df(db, sub, ann)
    print("nes_threshold=%.1f -> shape %s" % (nes, df.shape))
    if df.shape[0]:
        print("columns:", list(df.columns))
        ncol = [c for c in df.columns if "NES" in c][0]
        print(df.sort_values(by=ncol, ascending=False).head(8).to_string())
        break
