import re
s=open(r'D:\桌面\single cell\analysis_output\PMC13453878_europepmc.xml',encoding='utf-8').read()
for m in re.finditer(r'<graphic[^>]+xlink:href="([^"]+)"',s): print(m.group(1))
for m in re.finditer(r'<supplementary-material[^>]*>.*?</supplementary-material>',s,re.S):
 print(m.group(0)[:1000])
