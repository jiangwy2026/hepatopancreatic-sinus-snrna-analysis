import re
t=open(r'D:\桌面\single cell\analysis_output\PMC13453878_text.txt',encoding='utf-8').read()
for term in ['Data availability','Availability of data','GEO','GSE','PRJ','CNGB','SRA','accession','Supplementary Table S4','Table S4']:
 print('\n###',term)
 for m in list(re.finditer(re.escape(term),t,re.I))[:6]:
  print(t[max(0,m.start()-500):m.start()+1200])
