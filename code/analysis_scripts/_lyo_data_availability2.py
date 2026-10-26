import re
t=open(r'D:\桌面\single cell\analysis_output\PMC13453878_text.txt',encoding='utf-8').read()
out=[]
for term in ['Data availability','Availability of data','GEO','GSE','PRJ','CNGB','SRA','accession','Supplementary Table S4','Table S4']:
 out.append('\n### '+term)
 for m in list(re.finditer(re.escape(term),t,re.I))[:8]:
  out.append(t[max(0,m.start()-500):m.start()+1500])
open(r'D:\桌面\single cell\analysis_output\lyo_data_availability.txt','w',encoding='utf-8').write('\n'.join(out))
print('ok')
