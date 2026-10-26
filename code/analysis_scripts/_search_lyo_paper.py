import urllib.request, urllib.parse
q='Single-nucleus RNA sequencing illuminates a functional analog of lymphoid organ in crustacean'
url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term='+urllib.parse.quote(q)+'&retmode=json'
try:
 print(urllib.request.urlopen(url,timeout=30).read().decode()[:2000])
except Exception as e: print('ERR',repr(e))
