p = r"D:\桌面\single cell\analysis_output\scripts\89_velocity_scvelo.py"
src = open(p, encoding="utf8").read()
old = """scv.pp.filter_and_normalize(A, min_shared_counts=20)
sc.pp.highly_variable_genes(A, n_top_genes=2000)
print("HVGs:", int(A.var["highly_variable"].sum()))
scv.pp.moments(A, n_pcs=30, n_neighbors=30)"""
new = """scv.pp.filter_and_normalize(A, min_shared_counts=20)
print("genes kept:", A.n_vars)
# guard against non-finite values introduced upstream (SoupX correction can leave them)
for _k in (None, "spliced", "unspliced", "ambiguous"):
    _M = A.X if _k is None else A.layers.get(_k)
    if _M is None:
        continue
    if sp.issparse(_M):
        _bad = ~np.isfinite(_M.data)
        if _bad.any():
            _M.data[_bad] = 0.0
            print("  cleaned %d non-finite values in %s" % (int(_bad.sum()), _k))
    else:
        _bad = ~np.isfinite(_M)
        if _bad.any():
            _M[_bad] = 0.0
            print("  cleaned %d non-finite values in %s" % (int(_bad.sum()), _k))
scv.pp.moments(A, n_pcs=30, n_neighbors=30)"""
assert old in src
src = src.replace(old, new, 1)
open(p, "w", encoding="utf8").write(src)
print("patched")
