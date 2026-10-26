import urllib.request
url='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13453878/fullTextXML'
try:
 s=urllib.request.urlopen(url,timeout=60).read().decode('utf-8','replace'); print(len(s)); open(r'D:\桌面\single cell\analysis_output\PMC13453878_europepmc.xml','w',encoding='utf-8').write(s); print(s[:500])
except Exception as e: print('ERR',repr(e))
