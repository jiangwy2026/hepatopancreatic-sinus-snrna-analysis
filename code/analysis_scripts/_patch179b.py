p = r"D:\桌面\single cell\analysis_output\scripts\179_composition_stacked_preview.py"
s = open(p, encoding="utf-8").read()
s = s.replace('bbox_to_anchor=(0.50, 0.165)', 'bbox_to_anchor=(0.50, 0.205)')
old = '''    fig.text(0.50, 0.065,
             "PH 2.30×; MH2 2.17×; PM-like 0.27×   (treated / control)",
             ha="center", va="center", fontsize=5.2, color="#333333")
    fig.text(0.50, 0.018,
             "χ² = %.1f; P = 1.1 × 10$^{-135}$" % chi["chi2"],
             ha="center", va="center", fontsize=5.1, color="#333333")'''
s = s.replace(old, '')
open(p,"w",encoding="utf-8").write(s)
print("patched")
