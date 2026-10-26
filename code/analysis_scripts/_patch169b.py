p = r"D:\桌面\single cell\analysis_output\scripts\169_single_gene_tsne_10groups.py"
src = open(p, encoding="utf8").read().replace("\r\n", "\n")
old = '''        ax.set_title("%s [%s]\\n%.1f%% cells" % (sym, g, 100 * (v > 0).mean()), fontsize=4.9, pad=2)'''
new = '''        in_g = v[ct == g]
        ax.set_title("%s [%s]\\n%.0f%% of group / %.1f%% all" % (sym, g, 100 * (in_g > 0).mean(), 100 * (v > 0).mean()),
                     fontsize=4.7, pad=2)'''
assert old in src; src = src.replace(old, new, 1)
open(p, "w", encoding="utf8", newline="\n").write(src); print("patched")
