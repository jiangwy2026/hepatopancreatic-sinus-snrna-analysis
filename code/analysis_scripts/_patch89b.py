p = r"D:\桌面\single cell\analysis_output\scripts\89_velocity_scvelo.py"
src = open(p, encoding="utf8").read()
old = "scv.pp.filter_and_normalize(A, min_shared_counts=20, n_top_genes=2000)\nscv.pp.moments(A, n_pcs=30, n_neighbors=30)"
new = """scv.pp.filter_and_normalize(A, min_shared_counts=20)
sc.pp.highly_variable_genes(A, n_top_genes=2000)
print("HVGs:", int(A.var["highly_variable"].sum()))
scv.pp.moments(A, n_pcs=30, n_neighbors=30)"""
assert old in src
src = src.replace(old, new, 1)
open(p, "w", encoding="utf8").write(src)
print("patched")
