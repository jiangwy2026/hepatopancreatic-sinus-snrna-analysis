import re, html
s=open(r'D:\桌面\single cell\analysis_output\PMC13453878.xml',encoding='utf-8').read()
# crude text
t=re.sub(r'<[^>]+>',' ',s); t=html.unescape(t); t=re.sub(r'\s+',' ',t)
open(r'D:\桌面\single cell\analysis_output\PMC13453878_text.txt','w',encoding='utf-8').write(t)
for m in re.finditer(r'.{0,250}(?:cluster|Cluster|subcluster|cell type|cell-type).{0,350}',t):
 print('\n---\n',m.group(0)[:900])
