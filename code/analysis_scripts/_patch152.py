p = r"D:\桌面\single cell\analysis_output\scripts\152_bubble_split.py"
src = open(p, encoding="utf8").read().replace("\r\n", "\n")
old = '''    cb = fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.02)'''
new = '''    cb = fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.22)'''
assert old in src; src = src.replace(old, new, 1)
old2 = '''    for k, v in enumerate((5, 20, 50)):
        ax.scatter(len(cols) - 0.15, 0.6 + k * 1.15, s=10 + v * 3.6, c="none", edgecolors="#666666",
                   linewidths=0.4, clip_on=False)
        ax.text(len(cols) + 0.08, 0.6 + k * 1.15, "%d genes" % v, fontsize=5.6, va="center", clip_on=False)'''
new2 = '''    ax.text(len(cols) + 0.12, -0.30, "genes", fontsize=5.8, va="center", clip_on=False)
    for k, v in enumerate((5, 20, 50)):
        ax.scatter(len(cols) + 0.10, 0.35 + k * 1.05, s=10 + v * 3.6, c="none", edgecolors="#666666",
                   linewidths=0.4, clip_on=False)
        ax.text(len(cols) + 0.42, 0.35 + k * 1.05, "%d" % v, fontsize=5.8, va="center", clip_on=False)'''
assert old2 in src; src = src.replace(old2, new2, 1)
open(p, "w", encoding="utf8", newline="\n").write(src); print("patched")
