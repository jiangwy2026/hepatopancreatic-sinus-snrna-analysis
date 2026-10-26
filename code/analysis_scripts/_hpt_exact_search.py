import urllib.request, urllib.parse, json
qs={
 'exact_hpt_shrimp':'TITLE_ABS:("hematopoietic tissue" OR "haematopoietic tissue" OR "hemopoietic tissue") AND TITLE_ABS:(shrimp OR Penaeus OR Litopenaeus OR prawn)',
 'exact_hpt_vannamei':'TITLE_ABS:("hematopoietic tissue" OR "haematopoietic tissue" OR "hemopoietic tissue") AND TITLE_ABS:("Litopenaeus vannamei" OR "Penaeus vannamei")',
 'exact_HPT_shrimp':'TITLE_ABS:(HPT) AND TITLE_ABS:(shrimp OR Penaeus OR Litopenaeus)'
}
for k,q in qs.items():
 url='https://www.ebi.ac.uk/europepmc/webservices/rest/search?query='+urllib.parse.quote(q)+'&format=json&pageSize=30&resultType=core'
 try:
  d=json.loads(urllib.request.urlopen(url,timeout=60).read().decode('utf-8','replace'))
  print('\n###',k,'hits',d.get('hitCount'))
  for r in d.get('resultList',{}).get('result',[])[:20]:
   print(r.get('pubYear'),r.get('title'),'| PMID',r.get('pmid'),'| DOI',r.get('doi'))
 except Exception as e: print('ERR',k,repr(e))
