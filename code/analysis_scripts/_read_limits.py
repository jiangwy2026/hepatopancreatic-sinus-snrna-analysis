import io, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
p = r"D:\桌面\single cell\文章所需\cell report\_extract_MANUSCRIPTS_0702.txt"
paras = [x for x in open(p, encoding="utf8").read().split("\n") if x.strip()]
for i in (42, 43, 184, 185, 186, 187, 188):
    print("――― [%d] ―――" % i); print(paras[i]); print()
