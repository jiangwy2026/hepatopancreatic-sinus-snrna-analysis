# -*- coding: utf-8 -*-
"""185_scenic_prune_annotation_test.py - test whether the v9 motif annotation matches the local ranking DB."""
import os, sys, time, pickle
import numpy as np, pandas as pd
if not hasattr(np, "object"): np.object = object
import dask
dask.config.set({"dataframe.query-planning": False})
OUT=r"D:\桌面\single cell\analysis_output"; WORK=r"E:\workspace\scenic_work"; DB=os.path.join(OUT,"scenic_db")

def main():
    expr=pd.read_csv(os.path.join(WORK,"expr.csv.gz"),index_col=0)
    adj=pd.read_csv(os.path.join(OUT,"scenic_adjacencies.csv")); adj.columns=["TF","target","weight"]
    adj=adj[adj.TF.isin(expr.columns)&adj.target.isin(expr.columns)].reset_index(drop=True)
    from pyscenic.utils import modules_from_adjacencies
    mods=list(modules_from_adjacencies(adj,expr,rho_mask_dropouts=False))
    print("modules",len(mods),flush=True)
    from ctxcore.rnkdb import FeatherRankingDatabase as RDB
    from pyscenic.prune import prune2df,df2regulons
    rdb=RDB(os.path.join(DB,"dm6_v10_clust.genes_vs_motifs.rankings.feather"),"dm6")
    for tag,tbl in [("v9","motifs-v9-nr.flybase.tbl"),("v10","motifs-v10nr_clust-nr.flybase-m0.001-o0.0.tbl")]:
        t0=time.time()
        try:
            df=prune2df([rdb],mods,os.path.join(DB,tbl),num_workers=5)
            print(tag,"prune shape",df.shape,"time",round(time.time()-t0,1),flush=True)
            if not df.empty:
                df.to_csv(os.path.join(OUT,"scenic_motif_enrichment_%s.csv"%tag))
                regs=df2regulons(df)
                pickle.dump(regs,open(os.path.join(OUT,"scenic_regulons_%s.pkl"%tag),"wb"))
                print(tag,"regulons",len(regs),[r.name for r in regs[:20]],flush=True)
        except Exception as e:
            print(tag,"failed",repr(e)[:300],flush=True)
if __name__=="__main__":
    log=open(r"D:\桌面\single cell\logs\scenic_annotation_test.log","w",encoding="utf8",buffering=1)
    sys.stdout=log; sys.stderr=log; main(); log.close()
