p = r"D:\桌面\single cell\analysis_output\scripts\127_syncytium_tests.py"
src = open(p, encoding="utf8").read().replace("\r\n", "\n")
old = '''gp = {g: i for i, g in enumerate(raw.var_names.astype(str))}
def panel_scores(idx):'''
new = '''gp = {g: i for i, g in enumerate(raw.var_names.astype(str))}
# 符号 -> 取总 counts 最高的那个基因 ID（避免同名基因取错）
_sym = raw.var["sym"].astype(str).values
_tot_gene = np.asarray(X.sum(0)).ravel()
_best = {}
for i, s in enumerate(_sym):
    if s not in _best or _tot_gene[i] > _tot_gene[_best[s]]:
        _best[s] = i
def _ids(gl):
    return [gp[str(raw.var_names[_best[x]])] for x in gl if x in _best]
def panel_scores(idx):'''
assert old in src; src = src.replace(old, new, 1)
src = src.replace('''    for nm, gl in (("S", S_GENES), ("G2M", G2M_GENES)):
        gi = [gp[x] for x in gl if x in gp]''', '''    for nm, gl in (("S", S_GENES), ("G2M", G2M_GENES)):
        gi = _ids(gl)''')
assert "_ids(gl)" in src
open(p, "w", encoding="utf8", newline="\n").write(src); print("patched")
