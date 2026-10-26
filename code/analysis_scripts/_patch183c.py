p = r"D:\桌面\single cell\analysis_output\scripts\183_y_branch_diffmap_elife_preview.py"
s = open(p, encoding="utf-8").read()
old1 = '''    # Dashed lineage ellipses
    for g in ["PH", "GH", "MH2"]:
        mu, (w, h), ang = ell(D[bct == g], 1.45)
        ax.add_patch(Ellipse(mu, w, h, angle=ang, fill=False, ls=(0, (4, 2.2)),
                             lw=0.95, ec=COL[g], zorder=6))
    # central TS ring
    mu, (w, h), ang = ell(D[bct == "TS"], 1.1)
    ax.add_patch(Ellipse(mu, w, h, angle=ang, fill=False, ls=(0, (2, 2)),
                         lw=0.7, ec="#555555", zorder=5))
'''
s = s.replace(old1, '')
open(p,"w",encoding="utf-8").write(s)
print("patched")
