import re
from pathlib import Path
from docx import Document
p=Path(r"D:\桌面\single cell\修订交付_2026-10-26_v2\manuscripts_1026_revised_final_2026-10-26_final.docx")
d=Document(p)
paras=d.paragraphs
ref_start=next(i for i,x in enumerate(paras) if x.text.strip().startswith('References:'))
refs=[x.text for x in paras if x.style.name=='EndNote Bibliography']
first=[]; occ=[]
for i,x in enumerate(paras[:ref_start]):
    t=x.text
    for m in re.finditer(r'\(([0-9][0-9,\s–-]*)\)',t):
        s=m.group(1)
        if not re.fullmatch(r'[0-9,\s–-]+',s): continue
        nums=[int(v) for v in re.findall(r'\d+',s)]
        for n in nums:
            if 1<=n<=len(refs): occ.append((i,n))
            if n not in first: first.append(n)
print('refs',len(refs),'first_count',len(first),'missing',sorted(set(range(1,len(refs)+1))-set(first)))
print('first order',first)
print('out_of_order',[(first[i],first[i+1]) for i in range(len(first)-1) if first[i+1]<first[i]])
print('max',max(first),'occurrences',len(occ))
print('last refs',[(i,refs[i-1][:100]) for i in range(max(1,len(refs)-10),len(refs)+1)])
print('size',p.stat().st_size)
