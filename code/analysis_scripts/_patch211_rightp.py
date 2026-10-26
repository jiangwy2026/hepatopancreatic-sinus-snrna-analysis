p=r"D:\桌面\single cell\analysis_output\scripts\211_alluvial_supported_only.py"
s=open(p,encoding="utf-8").read()
# build right_p mapping before drawing nodes
s=s.replace(''' for o in ORGS:
  sub=S[S.organ==o].copy(); sub["label"]=sub.best_group.astype(str).map(CELL[o]).fillna(sub.best_group.astype(str)+" (unknown)")
  labs=list(dict.fromkeys(sub["label"].tolist())); labels[o]=labs
  ys=np.linspace(mid_y[o]+0.07,mid_y[o]-0.07,len(labs)) if len(labs)>1 else [mid_y[o]]
  for l,y in zip(labs,ys): right.append((o,l,y))
 rpos={(o,l):y for o,l,y in right}''',''' right_p={}
 for o in ORGS:
  sub=S[S.organ==o].copy(); sub["label"]=sub.best_group.astype(str).map(CELL[o]).fillna(sub.best_group.astype(str)+" (unknown)")
  labs=list(dict.fromkeys(sub["label"].tolist())); labels[o]=labs
  ys=np.linspace(mid_y[o]+0.07,mid_y[o]-0.07,len(labs)) if len(labs)>1 else [mid_y[o]]
  for l,y in zip(labs,ys):
   right.append((o,l,y))
   right_p[(o,l)]=float(sub[sub.label==l].p_perm.min())
 rpos={(o,l):y for o,l,y in right}''')
# right node p label
s=s.replace('''  ax.scatter([0.51],[y],s=16,c=COL[o],edgecolors="white",linewidths=0.4,zorder=8); ax.text(0.525,y,l,ha="left",va="center",fontsize=4.4,color=COL[o])''','''  ax.scatter([0.51],[y],s=16,c=COL[o],edgecolors="white",linewidths=0.4,zorder=8); ax.text(0.525,y,l,ha="left",va="center",fontsize=4.4,color=COL[o])
  ax.text(0.84,y,"p=%.3g"%right_p[(o,l)],ha="left",va="center",fontsize=4.0,color=COL[o])''')
# remove left p label block
old='''   # p label at left, stacked by link
   py=left_y[p]+(k-(n-1)/2)*0.022
   ax.text(0.115,py,"p=%.3g"%r.p_perm,ha="left",va="center",fontsize=4.0,color=COL[r.organ],
           bbox=dict(boxstyle="round,pad=0.10",fc="white",ec="none",alpha=0.72),zorder=9)'''
s=s.replace(old,'')
s=s.replace('fig=plt.figure(figsize=(5.4,4.0),facecolor="white")','fig=plt.figure(figsize=(5.7,4.0),facecolor="white")')
s=s.replace('Numbers: permutation p. All significant/trend links plus top-2 organ matches are shown; marker-level candidate sharing, not migration.','p values are shown at the right. Marker-level candidate sharing, not migration.')
s=s.replace('"Fig6_alluvial_supported_only_preview.png"','"Fig6_alluvial_supported_only_rightp_preview.png"')
open(p,"w",encoding="utf-8").write(s)
print("patched")
