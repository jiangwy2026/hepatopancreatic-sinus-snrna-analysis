p = r"D:\桌面\single cell\analysis_output\scripts\165_final_annotation_figs.py"
src = open(p, encoding="utf8").read().replace("\r\n", "\n")
src = src.replace("def draw(n_tab, a_tab, cols, title, fname, figsize, scale, legend, xfs=6.0, rot=45):",
                  "def draw(n_tab, a_tab, cols, title, fname, figsize, scale, legend, xfs=6.0, rot=45, lx=0.45):")
src = src.replace('    ax.text(len(cols) + 0.12, -0.30, "genes", fontsize=5.8, va="center", clip_on=False)',
                  '    ax.text(len(cols) + lx + 0.02, -0.30, "genes", fontsize=5.8, va="center", clip_on=False)')
src = src.replace('        ax.scatter(len(cols) + 0.10, 0.35 + k * 1.05, s=14 + v * scale, c="none", edgecolors="#666666",',
                  '        ax.scatter(len(cols) + lx, 0.35 + k * 1.05, s=14 + v * scale, c="none", edgecolors="#666666",')
src = src.replace('        ax.text(len(cols) + 0.42, 0.35 + k * 1.05, "%d" % v, fontsize=5.8, va="center", clip_on=False)',
                  '        ax.text(len(cols) + lx + 0.32, 0.35 + k * 1.05, "%d" % v, fontsize=5.8, va="center", clip_on=False)')
src = src.replace('''         (6.9, 3.4), 30, (1, 2, 4), xfs=5.8)''', '''         (6.9, 3.4), 30, (1, 2, 4), xfs=5.8, lx=1.15)''')
open(p, "w", encoding="utf8", newline="\n").write(src); print("patched")
