p = r"D:\桌面\single cell\analysis_output\scripts\183_y_branch_diffmap_elife_preview.py"
s = open(p, encoding="utf-8").read()
s = s.replace('ax.add_patch(Ellipse(mu, w, h, angle=ang, fill=False, ls=(0, (4, 2.2)),\n                             lw=1.0, ec=COL[g], zorder=6))', 'ax.add_patch(Ellipse(mu, w, h, angle=ang, fill=False, ls=(0, (4, 2.2)),\n                             lw=0.95, ec=COL[g], zorder=6))')
s = s.replace('mu, (w, h), ang = ell(D[bct == g], 1.65)', 'mu, (w, h), ang = ell(D[bct == g], 1.45)')
s = s.replace('mu, (w, h), ang = ell(D[bct == "TS"], 1.3)', 'mu, (w, h), ang = ell(D[bct == "TS"], 1.1)')
s = s.replace('    ax.set_title("Hemocyte lineage branches", fontsize=7.0, pad=18)\n', '')
s = s.replace('bbox_to_anchor=(0.5, 1.14)', 'bbox_to_anchor=(0.5, 1.10)')
open(p,"w",encoding="utf-8").write(s)
print("patched")
