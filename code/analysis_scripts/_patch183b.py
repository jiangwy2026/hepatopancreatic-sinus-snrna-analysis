p = r"D:\桌面\single cell\analysis_output\scripts\183_y_branch_diffmap_elife_preview.py"
s = open(p, encoding="utf-8").read()
# insert smooth helper after ell function
needle = '''def main():
    A = ad.read_h5ad(os.path.join(OUT, "adata_hemocyte_traj_geodesic.h5ad"))'''
helper = '''def bezier(p0, p2, bend=0.15, n=160):
    p0 = np.asarray(p0, float); p2 = np.asarray(p2, float)
    mid = (p0 + p2) / 2.0
    v = p2 - p0
    perp = np.array([-v[1], v[0]])
    p1 = mid + perp * bend
    t = np.linspace(0, 1, n)[:, None]
    return (1 - t) ** 2 * p0 + 2 * (1 - t) * t * p1 + t ** 2 * p2


def draw_arrowhead(ax, curve, color, lw=1.8):
    p_end = curve[-1]; p_prev = curve[-8]
    ax.add_patch(FancyArrowPatch(p_prev, p_end, arrowstyle="-|>,head_length=0.55,head_width=0.26",
                                 mutation_scale=11, lw=lw, color=color, alpha=0.98, zorder=9,
                                 shrinkA=0, shrinkB=0,
                                 path_effects=[pe.withStroke(linewidth=2.8, foreground="white")]))


def main():
    A = ad.read_h5ad(os.path.join(OUT, "adata_hemocyte_traj_geodesic.h5ad"))'''
s = s.replace(needle, helper)
old = '''    # subtle direction arrows
    arrows = [("PH", "TS", 0.05, "#333333", 1.0),
              ("TS", "GH", 0.10, COL["GH"], 1.0),
              ("TS", "MH2", -0.10, COL["MH2"], 1.0)]
    for a, b, rad, c, lw in arrows:
        ax.add_patch(FancyArrowPatch(cent[a], cent[b], connectionstyle="arc3,rad=%s" % rad,
                                     arrowstyle="-|>,head_length=0.45,head_width=0.22",
                                     mutation_scale=9, lw=lw, color=c, alpha=0.9,
                                     zorder=8, shrinkA=2, shrinkB=2,
                                     path_effects=[pe.withStroke(linewidth=2.0, foreground="white")]))'''
new = '''    # smooth shared trunk and branches
    trunk = bezier(cent["PH"], cent["TS"], bend=0.06)
    br_gh = bezier(cent["TS"], cent["GH"], bend=0.20)
    br_mh2 = bezier(cent["TS"], cent["MH2"], bend=-0.18)
    for curve, color, lw in [(trunk, "#333333", 2.0),
                             (br_gh, COL["GH"], 2.0),
                             (br_mh2, COL["MH2"], 2.0)]:
        ax.plot(curve[:, 0], curve[:, 1], lw=lw, color=color, zorder=8,
                solid_capstyle="round", path_effects=[pe.withStroke(linewidth=3.2, foreground="white")])
        draw_arrowhead(ax, curve, color, lw=lw)'''
s = s.replace(old, new)
open(p,"w",encoding="utf-8").write(s)
print("patched")
