import urllib.request
url='https://cdn.ncbi.nlm.nih.gov/pmc/blobs/3c5d/13453878/67708eaaf5e2/gr1.jpg'
b=urllib.request.urlopen(url,timeout=60).read()
open(r'D:\桌面\single cell\analysis_output\lyo_gr1.jpg','wb').write(b)
print(len(b))
