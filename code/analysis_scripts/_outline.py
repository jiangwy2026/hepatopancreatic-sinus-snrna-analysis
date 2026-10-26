import io, sys, re
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
p = r"D:\桌面\single cell\文章所需\cell report\_extract_MANUSCRIPTS_0702.txt"
txt = open(p, encoding="utf8").read()
paras = [x for x in txt.split("\n") if x.strip()]
print("总段落:", len(paras), "总字符:", len(txt))
print("=" * 100)
print("【疑似标题 / 小标题 / 图注】")
for i, x in enumerate(paras):
    s = x.strip()
    if len(s) < 120 and (re.match(r"^(Figure|Fig\.|Table|Supplementary|Extended Data)", s) or
                         re.match(r"^[A-Z][A-Za-z0-9 ,'\-\(\)/]{2,70}$", s) or
                         re.match(r"^\d+\.\s", s) or s.endswith(":") and len(s) < 60):
        if len(s) > 3:
            print("  [%3d] %s" % (i, s))
