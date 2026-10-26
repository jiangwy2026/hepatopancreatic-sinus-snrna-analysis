# -*- coding: utf-8 -*-
r"""90_scenic_finish.py - SCENIC cisTarget pruning + AUCell (E:\workspace\scenic_env\python.exe)."""
import os, sys, io, time, pickle
import numpy as np, pandas as pd
if not hasattr(np, "object"):      # numpy>=1.24 removed the alias that pyscenic still uses
    np.object = object
import dask
dask.config.set({"dataframe.query-planning": False})   # legacy dask.dataframe API without dask-expr

OUT = r"D:\桌面\single cell\analysis_output"
WORK = r"E:\workspace\scenic_work"
DB = os.path.join(OUT, "scenic_db")

def main():
    expr = pd.read_csv(os.path.join(WORK, "expr.csv.gz"), index_col=0)
    print("expr matrix:", expr.shape, flush=True)
    adj = pd.read_csv(os.path.join(OUT, "scenic_adjacencies.csv"))
    adj.columns = ["TF", "target", "importance"]
    adj = adj[adj["TF"].isin(expr.columns) & adj["target"].isin(expr.columns)].reset_index(drop=True)
    print("adjacencies used:", adj.shape, flush=True)

    from pyscenic.utils import modules_from_adjacencies
    mods = list(modules_from_adjacencies(adj, expr, rho_mask_dropouts=False))
    print("modules:", len(mods), flush=True)

    from ctxcore.rnkdb import FeatherRankingDatabase as RDB
    from pyscenic.prune import prune2df, df2regulons
    rdb = RDB(os.path.join(DB, "dm6_v10_clust.genes_vs_motifs.rankings.feather"), "dm6")
    motif_ann = os.path.join(DB, "motifs-v10nr_clust-nr.flybase-m0.001-o0.0.tbl")
    t0 = time.time()
    df = prune2df([rdb], mods, motif_ann, num_workers=5)
    print("prune2df done in %.0fs, shape %s" % (time.time() - t0, df.shape), flush=True)
    df.to_csv(os.path.join(OUT, "scenic_motif_enrichment.csv"))
    regs = df2regulons(df)
    pickle.dump(regs, open(os.path.join(OUT, "scenic_regulons.pkl"), "wb"))
    names = sorted(r.name for r in regs)
    with open(os.path.join(OUT, "scenic_regulon_names.txt"), "w") as f:
        f.write("\n".join(names))
    print("regulons:", len(regs), flush=True)
    print("regulon names:", names[:40], flush=True)

    from pyscenic.aucell import aucell
    t0 = time.time()
    auc = aucell(expr, regs, num_workers=5)
    auc.to_csv(os.path.join(OUT, "scenic_aucell.csv"))
    print("aucell done in %.0fs, shape %s" % (time.time() - t0, auc.shape), flush=True)

    try:
        cm = pd.read_csv(os.path.join(WORK, "cellmeta.csv"), index_col=0)
        print("cellmeta columns:", list(cm.columns)[:12], flush=True)
        gcol = None
        for c in ("group", "Group", "celltype", "cluster", "leiden"):
            if c in cm.columns:
                gcol = c; break
        if gcol:
            common = [c for c in auc.index if c in cm.index]
            m = auc.loc[common].groupby(cm.loc[common, gcol].astype(str)).mean()
            m.T.to_csv(os.path.join(OUT, "scenic_auc_by_group.csv"))
            print("mean AUC by group written; groups:", list(m.index), flush=True)
        else:
            print("no group column found in cellmeta", flush=True)
    except Exception as e:
        print("group summary skipped:", repr(e)[:150], flush=True)
    print("DONE")

if __name__ == "__main__":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    main()
