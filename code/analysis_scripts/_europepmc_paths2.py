s=open(r'D:\桌面\single cell\analysis_output\PMC13453878_europepmc.xml',encoding='utf-8').read()
out=[]
for key in ['ga1.jpg','gr1.jpg','mmc1.pdf']:
 i=s.find(key); out.append('\n'+key+' '+str(i)+'\n'+s[max(0,i-600):i+1000])
open(r'D:\桌面\single cell\analysis_output\europepmc_paths.txt','w',encoding='utf-8').write('\n'.join(out))
print('ok')
