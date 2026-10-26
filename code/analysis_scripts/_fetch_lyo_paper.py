import urllib.request, urllib.parse, json
pmid='41248710'
for db in ['pubmed','pmc']:
 url=f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&db={db}&id={pmid}&retmode=json'
 try:
  s=urllib.request.urlopen(url,timeout=30).read().decode()
  print('\nDB',db,s[:3000])
 except Exception as e: print('ERR',db,repr(e))
url=f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id={pmid}&retmode=text&rettype=abstract'
try:
 print('\nABSTRACT\n',urllib.request.urlopen(url,timeout=30).read().decode()[:8000])
except Exception as e: print('ERR abstract',repr(e))
