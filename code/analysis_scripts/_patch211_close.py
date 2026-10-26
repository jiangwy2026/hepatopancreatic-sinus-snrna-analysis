p=r"D:\桌面\single cell\analysis_output\scripts\211_alluvial_supported_only.py"
s=open(p,encoding="utf-8").read()
old='''  ax.scatter([0.51],[y],s=16,c=COL[o],edgecolors="white",linewidths=0.4,zorder=8); ax.text(0.525,y,l,ha="left",va="center",fontsize=4.4,color=COL[o])
  ax.text(0.84,y,"p=%.3g"%right_p[(o,l)],ha="left",va="center",fontsize=4.0,color=COL[o])'''
new='''  ax.scatter([0.51],[y],s=16,c=COL[o],edgecolors="white",linewidths=0.4,zorder=8)
  ax.text(0.525,y,"%s   p=%.3g"%(l,right_p[(o,l)]),ha="left",va="center",fontsize=4.4,color=COL[o])'''
if old not in s: raise SystemExit('right p block not found')
s=s.replace(old,new)
s=s.replace('"Fig6_alluvial_supported_only_rightp_preview.png"','"Fig6_alluvial_supported_only_rightp_close_preview.png"')
open(p,"w",encoding="utf-8").write(s)
print("patched")
