p = r"D:\桌面\single cell\analysis_output\scripts\139_annotation_evidence.py"
src = open(p, encoding="utf8").read().replace("\r\n", "\n")
old = '''    def heat(ax, df, title, xlabels):
        M = df[xlabels].astype(float).values'''
new = '''    def heat(ax, df, title, cols, xlabels):
        M = df[cols].astype(float).values'''
assert old in src; src = src.replace(old, new, 1)
src = src.replace('im1 = heat(axb, E, "Published hemocyte atlas\\n(Yang et al. 2022, eLife)", XLAB)',
                  'im1 = heat(axb, E, "Published hemocyte atlas\\n(Yang et al. 2022, eLife)", CL, XLAB)')
src = src.replace('im2 = heat(axc, G, "Published hepatopancreas atlas\\n(GSE306674)", XLAB)',
                  'im2 = heat(axc, G, "Published hepatopancreas atlas\\n(GSE306674)", CL, XLAB)')
open(p, "w", encoding="utf8", newline="\n").write(src); print("patched")
