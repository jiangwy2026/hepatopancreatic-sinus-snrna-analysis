import urllib.request, gzip, os
url='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/789/085/GCF_003789085.1_ASM378908v1/GCF_003789085.1_ASM378908v1_feature_table.txt.gz'
out=r'D:\桌面\single cell\analysis_output\GCF_003789085.1_feature_table.txt.gz'
try:
 urllib.request.urlretrieve(url,out); print('downloaded',os.path.getsize(out))
except Exception as e: print('ERR',repr(e))
