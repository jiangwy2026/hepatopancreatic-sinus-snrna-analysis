p = r"D:\桌面\single cell\analysis_output\scripts\104_blood_dotplot_tsne.py"
src = open(p, encoding="utf8").read()
src = src.replace('    idx = {r["names"]: i for i, r in enumerate(H.var_names.astype(str))}',
                  '    idx = {g: i for i, g in enumerate(H.var_names.astype(str))}')
old = '''    sig = mk[(~mk.is_mt) & (mk.pvals_adj < 0.05) & (mk.logfoldchanges > 0.5)]'''
new = '''    import re
    HK = re.compile(r"^(Rp[LS]\\d|rpl-|rps-|RpLP|eef|eIF|Act\\d|Act5C|Tub|Hsp\\d|hsc\\d|HSC70|PRSS)", re.I)
    mk["hk"] = mk["sym"].str.match(HK).fillna(False)
    sig = mk[(~mk.is_mt) & (~mk.hk) & (mk.pvals_adj < 0.05) & (mk.logfoldchanges > 0.5)]'''
assert old in src
src = src.replace(old, new, 1)
open(p, "w", encoding="utf8").write(src)
print("patched")
