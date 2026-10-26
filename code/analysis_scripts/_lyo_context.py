import re
t=open(r'D:\桌面\single cell\analysis_output\PMC13453878_text.txt',encoding='utf-8').read()
terms=['16 types','ECM-producing','macrophage-like','lymphocyte precursor-like','capsular','SLC-rich','epithelial','progenitors','Table 1','Fig. 1','Figure 1']
out=[]
for term in terms:
 for m in re.finditer(re.escape(term),t):
  out.append('\n### '+term+'\n'+t[max(0,m.start()-800):m.start()+1600])
open(r'D:\桌面\single cell\analysis_output\PMC13453878_cluster_context.txt','w',encoding='utf-8').write('\n'.join(out))
print('written',len(out))
