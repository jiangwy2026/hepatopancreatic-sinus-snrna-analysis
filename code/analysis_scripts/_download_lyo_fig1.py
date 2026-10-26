import urllib.request
urls=['https://www.ncbi.nlm.nih.gov/pmc/articles/PMC13453878/bin/ga1.jpg','https://pmc.ncbi.nlm.nih.gov/articles/PMC13453878/bin/ga1.jpg']
for u in urls:
 try:
  b=urllib.request.urlopen(u,timeout=30).read(); open(r'D:\桌面\single cell\analysis_output\lyo_fig1.jpg','wb').write(b); print('ok',u,len(b)); break
 except Exception as e: print('err',u,repr(e))
