p = r"D:\桌面\single cell\analysis_output\scripts\133_pertype_tsne.py"
src = open(p, encoding="utf8").read().replace("\r\n", "\n")
old = '''    A.obs["orig_kept"] = [str(x) in set(map(str, ad.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"), backed="r").obs_names)) for x in A.obs_names]'''
new = '''    _sx = ad.read_h5ad(os.path.join(OUT, "adata_soupx.h5ad"), backed="r")
    _keptset = set(map(str, _sx.obs_names))
    A.obs["orig_kept"] = [str(x) in _keptset for x in A.obs_names]
    del _sx'''
assert old in src
open(p, "w", encoding="utf8", newline="\n").write(src.replace(old, new, 1)); print("patched")
