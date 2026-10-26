import urllib.request, urllib.parse, json
queries={
 'hpt_shrimp':'(Litopenaeus OR Penaeus OR shrimp) AND ("hematopoietic tissue" OR "haematopoietic tissue" OR "hemopoietic tissue" OR HPT)',
 'lyo_shrimp':'(Litopenaeus OR Penaeus OR shrimp) AND ("lymphoid organ") AND (hematopoiesis OR hemocyte OR haemocyte)',
 'hpt_crustacean':'(crustacean OR decapod) AND ("hematopoietic tissue" OR "haematopoietic tissue" OR "hemopoietic tissue" OR HPT)',
 'sc_shrimp_hem':'(Litopenaeus OR Penaeus OR shrimp) AND ("single-cell" OR "single-nucleus") AND (hematopoietic OR haemopoietic OR hemocyte OR haemocyte)'
}
out={}
for k,q in queries.items():
 url='https://www.ebi.ac.uk/europepmc/webservices/rest/search?query='+urllib.parse.quote(q)+'&format=json&pageSize=25&resultType=core'
 try:
  d=json.loads(urllib.request.urlopen(url,timeout=60).read().decode('utf-8','replace'))
  out[k]=d
  print('\n###',k,'hits',d.get('hitCount'))
  for r in d.get('resultList',{}).get('result',[])[:15]:
   print(r.get('pubYear'),r.get('title'),'|',r.get('journalTitle'),'| PMID',r.get('pmid'),'| DOI',r.get('doi'))
 except Exception as e: print('ERR',k,repr(e))
open(r'D:\桌面\single cell\analysis_output\hpt_literature_search.json','w',encoding='utf-8').write(json.dumps(out,ensure_ascii=False,indent=1))
