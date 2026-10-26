p = r"D:\桌面\single cell\analysis_output\scripts\100_traj_lines_final.py"
src = open(p, encoding="utf8").read()
old = """    for path in LIN:
        p = np.vstack([cent[g] for g in path]); c = dict(((a, b), col) for a, b, col in LIN)[(path[0], path[-1])]"""
new = """    for spec in LIN:
        path = spec[:-1]; c = spec[-1]
        p = np.vstack([cent[g] for g in path])"""
assert old in src
src = src.replace(old, new, 1)
open(p, "w", encoding="utf8").write(src)
print("patched")
