import io, sys, os
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
import pypdf, pdfplumber
root = r"D:\桌面\single cell\文章所需\cell report"
p = os.path.join(root, "CELL-REPORTS-S-26-05162.pdf")
r = pypdf.PdfReader(p)
print("页数:", len(r.pages))
t0 = (r.pages[0].extract_text() or "")
print("--- 第 1 页 ---"); print(t0[:1800])
with pdfplumber.open(p) as pdf:
    print("--- 图片/表格数（概览）---")
    for i, pg in enumerate(pdf.pages[:6]):
        print("  page %d: images %d, tables %d" % (i+1, len(pg.images), len(pg.find_tables())))
