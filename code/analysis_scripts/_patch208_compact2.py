p=r"D:\桌面\single cell\analysis_output\scripts\208_three_column_alluvial.py"
s=open(p,encoding="utf-8").read()
s=s.replace('fig.text(0.5,0.045,"Top 2 organ matches per LyO programme. Link width = z-score; marker-level candidate sharing, not migration.",ha="center",fontsize=4.4,color="#555555")','fig.text(0.5,0.075,"Top 2 organ matches per LyO programme. Link width = z-score; marker-level candidate sharing, not migration.",ha="center",fontsize=4.4,color="#555555")')
s=s.replace('loc="lower center",bbox_to_anchor=(0.5,0.015)','loc="lower center",bbox_to_anchor=(0.5,0.010)')
open(p,"w",encoding="utf-8").write(s)
print("patched")
