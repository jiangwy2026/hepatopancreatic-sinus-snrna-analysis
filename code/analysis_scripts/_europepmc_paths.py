s=open(r'D:\桌面\single cell\analysis_output\PMC13453878_europepmc.xml',encoding='utf-8').read()
for key in ['ga1.jpg','gr1.jpg','mmc1.pdf']:
 i=s.find(key); print('\n',key,i); print(s[max(0,i-500):i+800])
