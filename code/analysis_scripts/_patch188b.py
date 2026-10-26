p=r"D:\桌面\single cell\analysis_output\scripts\188_ccnd2_TGM_mechanism_preview.py"
s=open(p,encoding="utf-8").read()
s=s.replace('    ax.text(0.5, 0.96, "Candidate mechanism", ha="center", va="top", fontsize=6.8)\n','')
open(p,"w",encoding="utf-8").write(s)
print("patched")
