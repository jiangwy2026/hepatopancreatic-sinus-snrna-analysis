p = r"D:\桌面\single cell\analysis_output\scripts\129_0h_vs_6h_syncytium.py"
src = open(p, encoding="utf8").read().replace("\r\n", "\n")
old = '''def ids(gl): return [gp[str(var_ids[best[x]])] for x in gl if x in best]
def mean_expr(gl):
    gi = ids(gl)
    return np.asarray(ln[:, gi].mean(1)).ravel() if gi else np.zeros(len(sample))'''
new = '''def ids(gl): return [gp[str(var_ids[best[x]])] for x in gl if x in best]      # 输入基因名
def ids_by_id(gl): return [gp[x] for x in gl if x in gp]                       # 输入基因 ID
def mean_expr(gl):
    gi = ids(gl)
    return np.asarray(ln[:, gi].mean(1)).ravel() if gi else np.zeros(len(sample))
def mean_expr_id(gl):
    gi = ids_by_id(gl)
    return np.asarray(ln[:, gi].mean(1)).ravel() if gi else np.zeros(len(sample))'''
assert old in src; src = src.replace(old, new, 1)
src = src.replace('S["PH_sig"] = mean_expr(SIGPH)', 'S["PH_sig"] = mean_expr_id(SIGPH)')
src = src.replace('for k, v in PANEL.items(): S[k] = mean_expr(v)', 'for k, v in PANEL.items(): S[k] = mean_expr_id(v)')
open(p, "w", encoding="utf8", newline="\n").write(src); print("patched")
