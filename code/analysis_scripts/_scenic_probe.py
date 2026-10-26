import numpy as np
if not hasattr(np, "object"):
    np.object = object
try:
    import dask
    dask.config.set({"dataframe.query-planning": False})
    import dask.dataframe
    print("dask.dataframe OK")
except Exception as e:
    print("dask.dataframe FAIL:", repr(e)[:200])
try:
    from pyscenic.prune import prune2df, df2regulons
    print("pyscenic.prune OK")
except Exception as e:
    print("pyscenic.prune FAIL:", repr(e)[:300])
try:
    from pyscenic.aucell import aucell
    print("pyscenic.aucell OK")
except Exception as e:
    print("pyscenic.aucell FAIL:", repr(e)[:300])
try:
    from ctxcore.rnkdb import FeatherRankingDatabase as RDB
    print("ctxcore OK")
except Exception as e:
    print("ctxcore FAIL:", repr(e)[:200])
