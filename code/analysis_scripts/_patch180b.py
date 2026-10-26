p = r"D:\桌面\single cell\analysis_output\scripts\180_pseudotime_tsne_all_subtypes_preview.py"
s = open(p, encoding="utf-8").read()
# stronger arrows
s = s.replace('''    paths = [(cent["PH"], cent["TS"], -0.10, "#2B2B2B", 1.0),
             (cent["TS"], cent["GH"], 0.28, "#2166AC", 0.85),
             (cent["TS"], cent["MH2"], -0.28, "#8C510A", 0.85)]''',
'''    paths = [(cent["PH"], cent["TS"], 0.04, "#111111", 1.9),
             (cent["TS"], cent["GH"], 0.12, "#2166AC", 1.5),
             (cent["TS"], cent["MH2"], -0.12, "#8C510A", 1.5)]''')
s = s.replace('''                                     arrowstyle="-|>,head_length=0.45,head_width=0.20",
                                     mutation_scale=8, lw=lw, color=col, alpha=0.85,
                                     zorder=9, shrinkA=2, shrinkB=2,
                                     path_effects=[pe.withStroke(linewidth=1.8, foreground="white")]))''',
'''                                     arrowstyle="-|>,head_length=0.62,head_width=0.30",
                                     mutation_scale=13, lw=lw, color=col, alpha=0.98,
                                     zorder=12, shrinkA=3, shrinkB=3,
                                     path_effects=[pe.withStroke(linewidth=3.0, foreground="white")]))''')
# remove star and root annotation
old_root = '''    # root
    ph = np.where(ct == "PH")[0]
    root = int(ph[np.argmin(((U[ph] - U[ph].mean(0)) ** 2).sum(1))])
    ax.scatter(U[root, 0], U[root, 1], s=48, marker="*", c="#FF3B30",
               edgecolors="white", linewidths=0.55, zorder=11)
    ax.annotate("root", xy=(U[root, 0], U[root, 1]), xytext=(0, 11), textcoords="offset points",
                fontsize=5.3, color="#B2182B", ha="center",
                arrowprops=dict(arrowstyle="-", lw=0.5, color="#B2182B"), zorder=12)'''
s = s.replace(old_root, '    # PH root is indicated by the PH label; no single-cell star is shown.')
# label PH as root
s = s.replace('''    offs = {"PH": (-13, 9), "GH": (0, 0), "TS": (0, 0), "MH2": (0, 0),''',
                 '''    offs = {"PH": (-20, 9), "GH": (0, 0), "TS": (0, 0), "MH2": (0, 0),''')
s = s.replace('''        ax.annotate(g, xy=(c[0], c[1]), xytext=offs[g], textcoords="offset points",
                    fontsize=5.7, weight="bold", ha="center", va="center",''',
'''        lab = "PH (root)" if g == "PH" else g
        ax.annotate(lab, xy=(c[0], c[1]), xytext=offs[g], textcoords="offset points",
                    fontsize=5.7, weight="bold", ha="center", va="center",''')
open(p,"w",encoding="utf-8").write(s)
print("patched")
