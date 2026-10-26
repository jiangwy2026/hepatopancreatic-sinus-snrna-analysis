import io, sys, re
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
p = r"D:\桌面\single cell\文章所需\cell report\_extract_MANUSCRIPTS_0702.txt"
paras = [x for x in open(p, encoding="utf8").read().split("\n") if x.strip()]
print("=== 图注标题 ===")
for i, x in enumerate(paras):
    if re.match(r"^(Figure|Fig\.|Table)\s*\d", x.strip()) and len(x) < 400:
        print("  [%3d] %s" % (i, x.strip()))
print()
print("=== 补充图表提及（Supplementary Figure / Table）===")
seen = set()
for x in paras:
    for m in re.findall(r"(Supplementary (?:Figure|Table|Data)\s*\d+[A-Za-z]?)", x):
        if m not in seen:
            seen.add(m)
print("  ", ", ".join(sorted(seen)))
