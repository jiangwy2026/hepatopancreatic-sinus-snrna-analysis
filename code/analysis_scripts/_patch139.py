import os
p = r"D:\桌面\single cell\analysis_output\scripts\139_annotation_evidence.py"
src = open(p, encoding="utf8").read().replace("\r\n", "\n")
old = """    E.columns = [str(c) for c in E.columns]; G.columns = [str(c) for c in G.columns]"""
new = """    E.columns = [str(c) for c in E.columns]; G.columns = [str(c) for c in G.columns]
    EN = {"PH1(原血细胞1)": "PH1 (prohemocyte 1)", "PH2(原血细胞2)": "PH2 (prohemocyte 2)",
          "GH1(颗粒1)": "GH1 (granular 1)", "GH2(颗粒2)": "GH2 (granular 2)",
          "MH1(单核1)": "MH1 (monocytic 1)", "MH2(单核2·巨噬样)": "MH2 (monocytic 2, macrophage-like)"}
    E.index = [EN.get(str(i), str(i)) for i in E.index]
    XLAB = ["%s (C%s)" % (NAME[c], c) for c in CL]"""
assert old in src; src = src.replace(old, new, 1)
src = src.replace('def heat(ax, df, title, xlabels):', 'def heat(ax, df, title, xlabels):')
src = src.replace('im1 = heat(axb, E, "Published hemocyte atlas\\n(Yang et al. 2022, eLife)", CL)',
                  'im1 = heat(axb, E, "Published hemocyte atlas\\n(Yang et al. 2022, eLife)", XLAB)')
src = src.replace('im2 = heat(axc, G, "Published hepatopancreas atlas\\n(GSE306674)", CL)',
                  'im2 = heat(axc, G, "Published hepatopancreas atlas\\n(GSE306674)", XLAB)')
src = src.replace('columns = our vendor clusters 0-9', 'columns = our clusters')
src = src.replace('lab += "\\n→ %s (%0.1f)" % (eb[0].split("(")[0], eb[1])', 'lab += "\\n→ %s (%0.1f)" % (str(eb[0]).split(" (")[0], eb[1])')
open(p, "w", encoding="utf8", newline="\n").write(src); print("patched")
