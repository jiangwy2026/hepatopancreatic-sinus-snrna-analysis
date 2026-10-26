import re
t=open(r'D:\桌面\single cell\analysis_output\PMC13453878_text.txt',encoding='utf-8').read()
for term in ['Fig. 1','Figure 1','Identification of cell types','cell types in Oka','major cell types']:
 print('\n###',term)
 for m in list(re.finditer(re.escape(term),t))[:5]:
  print(t[max(0,m.start()-200):m.start()+2500])
