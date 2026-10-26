p = r"D:\桌面\single cell\analysis_output\scripts\176_sc_summary_new_names_preview.py"
s = open(p, encoding="utf-8").read()
old = '''    ax.set_xlabel("Cells (%)", fontsize=6.0)
    ax.set_title("Cell composition (%)", fontsize=7.0, pad=5)
    ax.text(0.02, 1.018, "χ² = %.1f; $P$ = %.1f × 10$^{-135}$" % (chi["chi2"], chi["p"] / 1e-135),
            transform=ax.transAxes, fontsize=5.3, va="bottom")'''
new = '''    ax.set_xlabel("Cells (%%)   (χ² = %.1f; $P$ = %.1f × 10$^{-135}$)" % (chi["chi2"], chi["p"] / 1e-135),
                  fontsize=5.5)
    ax.set_title("Cell composition (%)", fontsize=7.0, pad=5)'''
s = s.replace(old, new)
open(p,"w",encoding="utf-8").write(s)
print("patched composition")
