import urllib.request
url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=13453878&retmode=xml'
s=urllib.request.urlopen(url,timeout=60).read().decode('utf-8','replace')
open(r'D:\桌面\single cell\analysis_output\PMC13453878.xml','w',encoding='utf-8').write(s)
print('len',len(s))
