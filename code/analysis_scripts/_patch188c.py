p=r"D:\桌面\single cell\analysis_output\scripts\188_ccnd2_TGM_mechanism_preview.py"
s=open(p,encoding="utf-8").read()
# add negative volcano labels
old = '''    labels = {
        "LOC113823945": ("hemocyte TGM", (-1.3, 2.0), "#B2182B"),
        "ccnd2": ("ccnd2", (1.0, 1.0), "#B2182B"),
        "SSPO": ("SSPO", (-0.2, 1.4), "#B2182B"),
        "LOC113808957": ("Notch-like", (0.4, 2.8), "#B2182B"),
    }'''
new = '''    labels = {
        "LOC113823945": ("hemocyte TGM", (-1.3, 2.0), "#B2182B"),
        "ccnd2": ("ccnd2", (1.0, 1.0), "#B2182B"),
        "SSPO": ("SSPO", (-0.2, 1.4), "#B2182B"),
        "LOC113808957": ("Notch-like", (0.4, 2.8), "#B2182B"),
        "THBS4": ("THBS4", (-0.8, 3.2), "#2166AC"),
        "PEN-3": ("PEN-3", (-1.0, 2.6), "#2166AC"),
        "LOC113801825": ("crustin", (-0.8, 2.4), "#2166AC"),
        "LOC113810905": ("WAP", (-1.4, 2.4), "#2166AC"),
    }'''
s=s.replace(old,new)
# replace model function content at the end
old2 = '''    ax.text(0.50, 0.015, "Candidate model; TGM and ccnd2 are parallel programmes, not a proven direct interaction.",
            ha="center", va="bottom", fontsize=4.8, color="#555555")'''
new2 = '''    # negative / differentiation arm, relieved in treated PH
    box(ax, 0.70, 0.015, 0.20, 0.13,
        "Differentiation / restraining arm\\nTHBS4↓  PEN-3↓  crustin↓  WAP↓",
        "#E8F0F8", "#4575B4", fs=4.7, tc="#1F4E79")
    # blunt inhibitory connection to PH expansion
    ax.plot([0.80, 0.80], [0.145, 0.305], color="#4575B4", lw=0.7, ls=(0, (2, 2)), zorder=4)
    ax.plot([0.785, 0.815], [0.305, 0.305], color="#4575B4", lw=1.2, zorder=5)
    ax.text(0.815, 0.215, "relieved", fontsize=4.6, color="#4575B4", ha="left", va="center")'''
s=s.replace(old2,new2)
# add footer at figure level before save
s=s.replace('''    p = os.path.join(PREV, "04b_ccnd2_TGM_mechanism_preview.png")''','''    fig.text(0.50, 0.012,
             "Candidate model; the differentiation/restraining arm is suppressed, not a proven direct interaction.",
             ha="center", va="bottom", fontsize=4.8, color="#555555")
    p = os.path.join(PREV, "04b_ccnd2_TGM_mechanism_preview.png")''')
open(p,"w",encoding="utf-8").write(s)
print("patched")
