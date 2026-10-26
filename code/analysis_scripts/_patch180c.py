p = r"D:\桌面\single cell\analysis_output\scripts\180_pseudotime_tsne_all_subtypes_preview.py"
s = open(p, encoding="utf-8").read()
s = s.replace('mutation_scale=13, lw=lw, color=col, alpha=0.98,', 'mutation_scale=11, lw=lw, color=col, alpha=0.98,')
s = s.replace('''    offs = {"PH": (-20, 9), "GH": (0, 0), "TS": (0, 0), "MH2": (0, 0),
            "F": (0, -9), "R": (12, 0), "B": (8, 4), "M": (-12, -4),
            "MyoEpi": (10, -7), "PM-like": (0, 6)}''',
'''    offs = {"PH": (-20, 9), "GH": (-18, -11), "TS": (0, -15), "MH2": (0, 8),
            "F": (0, -9), "R": (12, 0), "B": (8, 4), "M": (-12, -4),
            "MyoEpi": (10, -7), "PM-like": (0, 6)}''')
open(p,"w",encoding="utf-8").write(s)
print("patched")
