p = r"D:\桌面\single cell\analysis_output\scripts\181_pseudotime_branch_tsne_elife_preview.py"
s = open(p, encoding="utf-8").read()
s = s.replace('ax.text(x, y, "PH → TS", fontsize=5.2, color=C_TRUNK, ha="center", va="center",', 'ax.text(x - 0.8, y + 0.45, "PH → TS", fontsize=5.2, color=C_TRUNK, ha="center", va="center",')
open(p,"w",encoding="utf-8").write(s)
print("patched")
