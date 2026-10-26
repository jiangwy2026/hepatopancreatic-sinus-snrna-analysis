p = r"D:\桌面\single cell\analysis_output\scripts\176_sc_summary_new_names_preview.py"
s = open(p, encoding="utf-8").read()
s = s.replace("\\n", "\n")
old = '''    for g in ORDER:
        c = cent[g]
        ax.text(c[0], c[1], g, fontsize=6.0, weight="bold", ha="center", va="center", color="white",
                path_effects=[pe.withStroke(linewidth=2.0, foreground="#333333")], zorder=10)'''
new = '''    offs = {"PH": (-12, 9), "GH": (0, 0), "TS": (0, 0), "MH2": (0, 0)}
    for g in ORDER:
        c = cent[g]
        ax.annotate(g, xy=(c[0], c[1]), xytext=offs[g], textcoords="offset points",
                    fontsize=6.0, weight="bold", ha="center", va="center", color="white",
                    path_effects=[pe.withStroke(linewidth=2.0, foreground="#333333")], zorder=10)'''
s = s.replace(old, new)
open(p, "w", encoding="utf-8").write(s)
print("patched")
