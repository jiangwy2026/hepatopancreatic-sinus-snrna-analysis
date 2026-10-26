import urllib.request, json
url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=bioproject&term=PRJNA1478458&retmode=json'
try:
 print(urllib.request.urlopen(url,timeout=20).read().decode()[:2000])
except Exception as e:
 print('ERR',repr(e))
