p = r"D:\桌面\single cell\analysis_output\scripts\89_velocity_scvelo.py"
src = open(p, encoding="utf8").read()
old = '        vpt.set_iroot(A.obs["is_root"].values.astype(bool))'
new = '        vpt.set_iroot("is_root")   # column name (bool array is not accepted by scvelo)'
assert old in src; src = src.replace(old, new, 1)
old2 = "    flow = np.asarray((Tr @ VG).todense()) / np.maximum(ncell_g[:, None], 1)"
new2 = "    flow = np.asarray((Tr @ VG @ Tr.T).todense()) / np.maximum(ncell_g[:, None], 1)"
assert old2 in src; src = src.replace(old2, new2, 1)
open(p, "w", encoding="utf8").write(src)
print("patched")
