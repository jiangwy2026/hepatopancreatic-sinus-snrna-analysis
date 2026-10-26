p = r"D:\桌面\single cell\analysis_output\scripts\179_composition_stacked_preview.py"
s = open(p, encoding="utf-8").read()
s = s.replace('fig = plt.figure(figsize=(4.8, 4.15), facecolor="white")', 'fig = plt.figure(figsize=(4.8, 3.55), facecolor="white")')
s = s.replace('ax = fig.add_axes([0.14, 0.40, 0.82, 0.42])', 'ax = fig.add_axes([0.14, 0.37, 0.82, 0.50])')
old = '''    handles = [Patch(facecolor=COL[g], edgecolor="white", label="%s (%s)" % (g, format(int(comp.loc[comp.group == g, "n_total"].iloc[0]), ",")))
               for g in order]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.50, 0.175),
               ncol=5, frameon=False, fontsize=5.1, handlelength=0.9,
               handletextpad=0.28, columnspacing=0.85, labelspacing=0.35)'''
new = '''    # Matplotlib fills legend columns first; interleave for row-wise reading.
    half = (len(order) + 1) // 2
    legend_order = []
    for i in range(half):
        legend_order.append(order[i])
        if i + half < len(order):
            legend_order.append(order[i + half])
    handles = [Patch(facecolor=COL[g], edgecolor="white", label="%s (%s)" % (g, format(int(comp.loc[comp.group == g, "n_total"].iloc[0]), ",")))
               for g in legend_order]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.50, 0.165),
               ncol=5, frameon=False, fontsize=5.1, handlelength=0.9,
               handletextpad=0.28, columnspacing=0.85, labelspacing=0.35)'''
s = s.replace(old, new)
s = s.replace('fig.text(0.50, 0.075,', 'fig.text(0.50, 0.065,')
s = s.replace('fig.text(0.50, 0.025,', 'fig.text(0.50, 0.018,')
open(p,"w",encoding="utf-8").write(s)
print("patched")
