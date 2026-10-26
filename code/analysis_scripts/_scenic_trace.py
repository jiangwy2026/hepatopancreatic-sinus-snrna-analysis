# -*- coding: utf-8 -*-
import os, numpy as np, pandas as pd
np.object = object
import dask; dask.config.set({"dataframe.query-planning": False})
from pyscenic import transform as TR
from pyscenic.utils import modules_from_adjacencies, load_motif_annotations
from ctxcore.rnkdb import FeatherRankingDatabase as RDB
OUT = r"D:\桌面\single cell\analysis_output"; DB = os.path.join(OUT, "scenic_db")
expr = pd.read_csv(r"E:\workspace\scenic_work\expr.csv.gz", index_col=0)
adj = pd.read_csv(os.path.join(OUT, "scenic_adjacencies.csv")); adj.columns = ["TF","target","importance"]
adj = adj[adj.TF.isin(expr.columns) & adj.target.isin(expr.columns)]
mods = list(modules_from_adjacencies(adj, expr, rho_mask_dropouts=False))
db = RDB(os.path.join(DB, "dm6_v10_clust.genes_vs_motifs.rankings.feather"), "dm6")
mod = mods[0]
print("module:", mod.name, "| n genes:", len(mod.genes), "| TF:", mod.transcription_factor)
df = db.load(mod)
print("db.load -> shape", df.shape, "| index name", df.index.name, "| first idx", list(df.index[:3]))
print("cols:", list(df.columns[:6]))
print("values[0,:6]:", df.values[0, :6], "dtype", df.values.dtype)
print("any NaN:", bool(pd.isnull(df.values).any()))
aucs = TR.calc_aucs(df, db.total_genes, np.ones(df.shape[1]), 0.05)
ness = (aucs - aucs.mean()) / aucs.std()
print("aucs min/median/max:", float(np.min(aucs)), float(np.median(aucs)), float(np.max(aucs)))
print("ness max:", float(np.max(ness)), "| n with NES>=3:", int((ness >= 3).sum()))
top = np.argsort(-ness)[:5]
print("top features:", [(str(df.index[i]), round(float(ness[i]), 2)) for i in top])
