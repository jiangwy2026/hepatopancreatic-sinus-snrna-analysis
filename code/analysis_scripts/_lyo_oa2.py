import urllib.request
for id in ['13453878','PMC13453878']:
 url='https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id='+id
 try:
  s=urllib.request.urlopen(url,timeout=30).read().decode('utf-8','replace'); print(id,s[:2000])
 except Exception as e: print(id,'ERR',repr(e))
