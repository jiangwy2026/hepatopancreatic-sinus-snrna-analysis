import docx, re, sys, os
sys.stdout.reconfigure(encoding="utf-8")
files = [
 r"D:\桌面\single cell\文章所需\cell report\MANUSCRIPTS_0702.docx",
 r"D:\桌面\single cell\文章所需\manuscripts_0617.docx",
 r"D:\桌面\single cell\文章所需\图注图标题.docx",
]
kw = ["niche","bone marrow","hematopoie","haematopoie","progenitor","sinusoid","stroma","support","microenvironment",
      "龛","骨髓","造血","前体","祖细胞","支持","微环境","血窦","基质"]
for p in files:
    print("="*80); print(p, os.path.exists(p))
    if not os.path.exists(p): continue
    d = docx.Document(p)
    paras = [x.text.strip() for x in d.paragraphs if x.text.strip()]
    print("paragraphs", len(paras), "chars", sum(map(len,paras)))
    hits = [i for i,t in enumerate(paras) if any(k.lower() in t.lower() for k in kw)]
    print("hit paragraphs", len(hits))
    for i in hits[:40]:
        print(f"\n[{i}] {paras[i][:600]}")
