import docx, sys
sys.stdout.reconfigure(encoding="utf-8")
p=r"D:\桌面\single cell\文章所需\cell report\MANUSCRIPTS_0702.docx"
d=docx.Document(p)
paras=[x.text.strip() for x in d.paragraphs if x.text.strip()]
for i in list(range(16,34))+list(range(140,158)):
    if i < len(paras):
        print(f"\n===== [{i}] =====")
        print(paras[i][:1800])
