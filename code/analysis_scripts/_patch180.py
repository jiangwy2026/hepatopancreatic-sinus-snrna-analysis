p = r"D:\桌面\single cell\analysis_output\scripts\180_pseudotime_tsne_all_subtypes_preview.py"
s = open(p, encoding="utf-8").read()
s = s.replace('fig = plt.figure(figsize=(3.95, 3.90), facecolor="white")', 'fig = plt.figure(figsize=(4.15, 3.90), facecolor="white")')
s = s.replace('ax = fig.add_axes([0.13, 0.105, 0.78, 0.78])', 'ax = fig.add_axes([0.115, 0.105, 0.74, 0.78])')
s = s.replace('    cax = ax.inset_axes([0.865, 0.56, 0.026, 0.24])', '    cax = fig.add_axes([0.875, 0.54, 0.020, 0.23])')
open(p,"w",encoding="utf-8").write(s)
print("patched")
