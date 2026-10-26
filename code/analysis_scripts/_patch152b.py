p = r"D:\桌面\single cell\analysis_output\scripts\152_bubble_split.py"
src = open(p, encoding="utf8").read().replace("\r\n", "\n")
# 面板独立尺寸标度
old = '''def draw(cols, title, fname, figsize, rot=45, xfs=6.2):
    n = N[cols].values; a = A[cols].values'''
new = '''def draw(cols, title, fname, figsize, rot=45, xfs=6.2, scale=3.6, legend=(5, 20, 50), lmax=None):
    n = N[cols].values; a = A[cols].values'''
assert old in src; src = src.replace(old, new, 1)
src = src.replace("ax.scatter(j, i, s=10 + n[i, j] * 3.6,", "ax.scatter(j, i, s=14 + n[i, j] * scale,")
src = src.replace("for k, v in enumerate((5, 20, 50)):", "for k, v in enumerate(legend):")
src = src.replace("ax.scatter(len(cols) + 0.10, 0.35 + k * 1.05, s=10 + v * 3.6, c=\"none\", edgecolors=\"#666666\",",
                  "ax.scatter(len(cols) + 0.10, 0.35 + k * 1.05, s=14 + v * scale, c=\"none\", edgecolors=\"#666666\",")
src = src.replace('''    draw(hps, "Hepatopancreas cell types (PMID 42004253)", "Fig_annotation_bubbles_hepatopancreas_preview.png",
         (6.6, 3.4), xfs=5.8)''',
'''    draw(hps, "Hepatopancreas cell types (PMID 42004253; bubble size scaled to 1-4 markers)",
         "Fig_annotation_bubbles_hepatopancreas_preview.png", (6.6, 3.4), xfs=5.8, scale=30, legend=(1, 2, 4))''')
open(p, "w", encoding="utf8", newline="\n").write(src); print("patched")
