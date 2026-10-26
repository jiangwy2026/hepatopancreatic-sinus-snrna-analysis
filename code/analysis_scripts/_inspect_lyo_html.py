s=open(r'D:\桌面\single cell\analysis_output\lyo_fig1.jpg',encoding='utf-8',errors='ignore').read()
for key in ['ga1.jpg','cdn.ncbi','image','href']:
 i=s.find(key); print('\n',key,i); print(s[max(0,i-300):i+800] if i>=0 else '')
