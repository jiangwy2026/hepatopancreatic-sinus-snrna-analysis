import urllib.request
url='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/789/085/GCF_003789085.1_ASM378908v1/'
try:
 s=urllib.request.urlopen(url,timeout=30).read().decode('utf-8','ignore')
 open(r'D:\桌面\single cell\analysis_output\GCF_ref_files.html','w',encoding='utf-8').write(s)
 print(s[:5000])
except Exception as e: print('ERR',repr(e))
