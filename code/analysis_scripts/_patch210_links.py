p=r"D:\桌面\single cell\analysis_output\scripts\210_alluvial_with_stats.py"
s=open(p,encoding="utf-8").read()
old='''  sub=P[(P.program==p)&(P.organ.isin(ORGS))].dropna(subset=["score"]).sort_values("score",ascending=False).head(2)
  for k,r in enumerate(sub.itertuples()):'''
new='''  allsub=P[(P.program==p)&(P.organ.isin(ORGS))].dropna(subset=["score"]).sort_values("score",ascending=False)
  sig=allsub[allsub.q_perm<0.10]
  if len(sig)>=2:
   sub=sig.copy()
  else:
   sub=pd.concat([sig,allsub.head(2)],ignore_index=True).drop_duplicates(subset=["organ"]).sort_values("score",ascending=False)
  nlinks=len(sub)
  for k,r in enumerate(sub.itertuples()):'''
if old not in s: raise SystemExit('selection block not found')
s=s.replace(old,new)
s=s.replace('yy=left_y[p]+(0.018 if k==0 else -0.018)','yy=left_y[p]+(k-(nlinks-1)/2)*0.018')
open(p,"w",encoding="utf-8").write(s)
print("patched")
