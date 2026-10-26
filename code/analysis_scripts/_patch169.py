p = r"D:\桌面\single cell\analysis_output\scripts\169_single_gene_tsne_10groups.py"
src = open(p, encoding="utf8").read().replace("\r\n", "\n")
old = '''    gs = GridSpec(nrow, ncol * 2, figure=fig, hspace=0.42, wspace=0.10,
                  left=0.02, right=0.985, top=0.93, bottom=0.03)'''
new = '''    gs = GridSpec(nrow, ncol, figure=fig, hspace=0.46, wspace=0.40,
                  left=0.02, right=0.975, top=0.93, bottom=0.03)'''
assert old in src; src = src.replace(old, new, 1)
old2 = '''        ax = fig.add_subplot(gs[r, c * 2]); cax = fig.add_subplot(gs[r, c * 2 + 1])'''
new2 = '''        ax = fig.add_subplot(gs[r, c])'''
assert old2 in src; src = src.replace(old2, new2, 1)
old3 = '''            cb = fig.colorbar(sm, cax=cax); cb.set_ticks([0, 1]); cb.ax.tick_params(labelsize=3.6, length=0.8)'''
new3 = '''            cax = ax.inset_axes([1.035, 0.30, 0.045, 0.40])
            cb = fig.colorbar(sm, cax=cax); cb.set_ticks([0, 1]); cb.ax.tick_params(labelsize=3.4, length=0.7, pad=0.8)'''
assert old3 in src; src = src.replace(old3, new3, 1)
open(p, "w", encoding="utf8", newline="\n").write(src); print("patched")
