# -*- coding: utf-8 -*-
import os, numpy as np, pandas as pd
from ctxcore.rnkdb import FeatherRankingDatabase as RDB
OUT = r"D:\桌面\single cell\analysis_output"; DB = os.path.join(OUT, "scenic_db")
expr = pd.read_csv(r"E:\workspace\scenic_work\expr.csv.gz", index_col=0)
for f in ("dm6_v10_clust.genes_vs_motifs.rankings.feather", "dm6_v10_clust.regions_vs_motifs.rankings.feather"):
    p = os.path.join(DB, f)
    try:
        db = RDB(p, "dm6")
        print("===", f)
        print("  total_genes:", db.total_genes, "| name:", db.name)
        try:
            print("  n motifs:", len(db.motifs), "| first:", list(db.motifs)[:5])
        except Exception as e:
            print("  motifs attr err:", repr(e)[:100])
        g = set(db.genes) if hasattr(db, "genes") else None
        if g:
            print("  genes in db:", len(g), "| overlap with expr:", len(g & set(expr.columns)))
            print("  example db genes:", sorted(g)[:8])
    except Exception as e:
        print("===", f, "FAILED:", repr(e)[:200])
for f in ("motifs-v10nr_clust-nr.flybase-m0.001-o0.0.tbl", "motifs-v9-nr.flybase.tbl"):
    p = os.path.join(DB, f)
    try:
        t = pd.read_csv(p, sep="\t", low_memory=False)
        print("===", f, t.shape)
        print("  columns:", list(t.columns))
        print("  head:", t.head(2).to_dict("records"))
    except Exception as e:
        print("===", f, "FAILED:", repr(e)[:200])
