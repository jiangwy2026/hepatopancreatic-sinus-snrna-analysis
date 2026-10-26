import re
s=open(r'D:\桌面\single cell\analysis_output\PMC13453878.xml',encoding='utf-8').read()
for m in re.finditer(r'<graphic[^>]+xlink:href="([^"]+)"',s):
 print(m.group(1))
