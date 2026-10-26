p = r"D:\桌面\single cell\analysis_output\scripts\124_dropped_bigcells.py"
src = open(p, encoding="utf8").read().replace("\r\n", "\n")
old = '''kept_idx = np.where(raw.obs["kept"].values)[0]
ct9 = raw.obs["vendor_cluster"].astype(str).map(BASE).replace({"NE": "NEU", "Neu": "NEU"}).values
ref = score(kept_idx); ref["ct"] = ct9[kept_idx]'''
new = '''kept_idx = np.where(raw.obs["kept"].values)[0]
lab = sx.obs["vendor_cluster"].astype(str).map(BASE).replace({"NE": "NEU", "Neu": "NEU"})
lab = lab.reindex([str(x) for x in raw.obs_names[kept_idx]]).values
ref = score(kept_idx); ref["ct"] = lab'''
assert old in src
open(p, "w", encoding="utf8", newline="\n").write(src.replace(old, new, 1))
print("patched")
