import re
t=open(r'D:\桌面\single cell\analysis_output\PMC13453878_text.txt',encoding='utf-8').read()
patterns=[r'Cluster\s+\d+[^.]{0,300}\.',r'cluster\s+\d+[^.]{0,300}\.',r'clusters?\s+\d+[^.]{0,300}\.',r'[^.]{0,200}cluster\s+\d+[^.]{0,250}\.']
seen=set(); out=[]
for pat in patterns:
 for m in re.finditer(pat,t):
  s=re.sub(r'\s+',' ',m.group(0)).strip()
  if len(s)>20 and s not in seen:
   seen.add(s); out.append(s)
open(r'D:\桌面\single cell\analysis_output\PMC13453878_cluster_sentences.txt','w',encoding='utf-8').write('\n'.join(out))
print('sentences',len(out))
