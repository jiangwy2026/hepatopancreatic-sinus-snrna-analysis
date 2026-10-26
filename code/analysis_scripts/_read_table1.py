import io, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
import pdfplumber, openpyxl, os
root = r"D:\桌面\single cell\文章所需\cell report"
print("=" * 90); print("TABLE 1 (table 1.pdf)"); print("=" * 90)
with pdfplumber.open(os.path.join(root, "table 1.pdf")) as pdf:
    for i, pg in enumerate(pdf.pages):
        t = pg.extract_text() or ""
        print("--- page %d ---" % (i + 1)); print(t)
        for tb in (pg.extract_tables() or []):
            for row in tb:
                print(" | ".join("" if c is None else str(c).replace("\n", " ") for c in row))
print(); print("=" * 90); print("KEY RESOURCES TABLE.xlsx"); print("=" * 90)
wb = openpyxl.load_workbook(os.path.join(root, "KEY RESOURCES TABLE.xlsx"), data_only=True)
for ws in wb.worksheets:
    print("--- sheet:", ws.title, ws.max_row, "rows ---")
    for r in ws.iter_rows(values_only=True):
        vals = ["" if v is None else str(v) for v in r]
        if any(v.strip() for v in vals):
            print(" | ".join(vals)[:300])
