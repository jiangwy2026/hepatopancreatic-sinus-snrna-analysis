p=r"D:\桌面\single cell\analysis_output\scripts\210_alluvial_with_stats.py"
s=open(p,encoding="utf-8").read()
old='''   # p label near mid-link
   yy=y[40]+(0.014 if r.organ=="Hepatic sinus" else -0.014 if r.organ=="Hemolymph" else 0.006)
   ax.text(0.185,yy,"p=%.3g"%r.p_perm,ha="center",va="center",fontsize=3.5,color=COL[r.organ],
           bbox=dict(boxstyle="round,pad=0.12",fc="white",ec="none",alpha=0.72),zorder=8)'''
new='''   # p label near the left node, stacked per programme to avoid mid-link collisions
   yy=left_y[p]+(0.018 if k==0 else -0.018)
   ax.text(0.115,yy,"p=%.3g"%r.p_perm,ha="left",va="center",fontsize=3.5,color=COL[r.organ],
           bbox=dict(boxstyle="round,pad=0.10",fc="white",ec="none",alpha=0.72),zorder=9)'''
if old not in s: raise SystemExit('p block not found')
s=s.replace(old,new)
s=s.replace('Line2D([0],[0],color="#999999",lw=0.9,ls=(0,(1,2)),alpha=0.35,label="unstable")]','Line2D([0],[0],color="#999999",lw=0.9,ls=(0,(1,2)),alpha=0.35,label="unstable"),\n                Line2D([0],[0],color="#777777",lw=0.8,alpha=0.25,label="stable n.s.")]')
open(p,"w",encoding="utf-8").write(s)
print("patched")
