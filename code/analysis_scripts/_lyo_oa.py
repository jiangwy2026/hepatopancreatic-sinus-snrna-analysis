import urllib.request
url='https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id=PMC13453878'
s=urllib.request.urlopen(url,timeout=30).read().decode('utf-8','replace')
open(r'D:\桌面\single cell\analysis_output\PMC13453878_oa.xml','w',encoding='utf-8').write(s)
print(s)
