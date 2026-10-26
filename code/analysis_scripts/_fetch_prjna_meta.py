import urllib.request, json, csv, io
url='https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1478458&result=read_run&fields=run_accession,sample_title,sample_alias,library_strategy,library_layout,fastq_ftp,fastq_bytes,read_count,base_count&format=tsv'
try:
 data=urllib.request.urlopen(url,timeout=60).read().decode('utf-8')
 open(r'D:\桌面\single cell\analysis_output\PRJNA1478458_ENA_metadata.tsv','w',encoding='utf-8').write(data)
 print(data[:4000])
 print('\nlines',len(data.splitlines()))
except Exception as e: print('ERR',repr(e))
