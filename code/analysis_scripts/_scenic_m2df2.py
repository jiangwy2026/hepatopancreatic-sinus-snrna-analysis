# -*- coding: utf-8 -*-
import os, numpy as np, pandas as pd
np.object = object
import dask; dask.config.set({"dataframe.query-planning": False})
from pyscenic import prune as P
from pyscenic.utils import load_motif_annotations, modules_from_adjacencies
from ctxcore.rnkdb import FeatherRankingDatabase as RDB
OUT = r"D:\桌面\single cell\analysis_output"; DB = os.path.join(OUT, "scenic_db")
ann_path = os.path.join(DB, "motifs-v10nr_clust-nr.flybase-m0.001-o0.0.tbl")
ann = load_motif_annotations(ann_path)
print("annotation frame:", ann.shape, "| index names:", ann.index.names, "| cols:", list(ann.columns)[:8])
print(ann.head(3).to_string()[:400])
expr = pd.read_csv(r"E:\workspace\scenic_work\expr.csv.gz", index_col=0)
adj = pd.read_csv(os.path.join(OUT, "scenic_adjacencies.csv")); adj.columns = ["TF","target","importance"]
adj = adj[adj.TF.isin(expr.columns) & adj.target.isin(expr.columns)]
mods = list(modules_from_adjacencies(adj, expr, rho_mask_dropouts=False))
db = RDB(os.path.join(DB, "dm6_v10_clust.genes_vs_motifs.rankings.feather"), "dm6")
sub = mods[:3]
df = P.modules2df(db, sub, ann)
print("modules2df raw ->", df.shape)
if df.shape[0]:
    print("index names:", df.index.names, "| columns:", list(df.columns))
    print(df.sort_values("NES", ascending=False).head(6).to_string()[:600])
