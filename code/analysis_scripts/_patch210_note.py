p=r"D:\桌面\single cell\analysis_output\scripts\210_alluvial_with_stats.py"
s=open(p,encoding="utf-8").read()
s=s.replace('fontsize=3.5,color=COL[r.organ]','fontsize=4.0,color=COL[r.organ]')
s=s.replace('"Numbers: permutation p. Top 2 organ matches per programme; marker-level candidate sharing, not migration."','"Numbers: permutation p. All significant/trend links plus top-2 organ matches are shown; marker-level candidate sharing, not migration."')
open(p,"w",encoding="utf-8").write(s)
print("patched")
