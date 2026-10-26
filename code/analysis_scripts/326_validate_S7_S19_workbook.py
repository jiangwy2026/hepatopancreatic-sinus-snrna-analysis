import openpyxl, json
from pathlib import Path
p=Path(r"D:\桌面\single cell\修订交付_2026-10-26_v2\Supplementary_Tables_S7-S19_2026-10-26.xlsx")
wb=openpyxl.load_workbook(p,read_only=True,data_only=True)
issues=[]
for ws in wb.worksheets:
    if ws.title=="README": continue
    # Read all values; validate duplicate/blank IDs where relevant.
    rows=ws.iter_rows(min_row=6,values_only=True)
    hdr=next(rows)
    data=list(rows)
    if not data: issues.append(f"{ws.title}: no data")
    if ws.title.startswith("S19_"):
        ids=[r[0] for r in data]
    elif "nucleus_id" in hdr:
        i=hdr.index("nucleus_id"); ids=[r[i] for r in data]
    elif "sample_id" in hdr:
        i=hdr.index("sample_id"); ids=[r[i] for r in data]
    elif "gene_id" in hdr:
        i=hdr.index("gene_id"); ids=[r[i] for r in data]
    else:
        ids=[]
    if ids and len(ids)!=len(set(ids)) and ws.title!="S19_Fig1":
        # marker/DEG IDs repeat across groups; only warn for nucleus/sample tables
        if ws.title in {"S7","S8","S9"}:
            issues.append(f"{ws.title}: duplicate primary IDs")
    print(ws.title,"rows",len(data),"cols",len(hdr))
wb.close()
# Cross-table count checks using source JSON.
base=Path(r"D:\桌面\single cell\修订交付_2026-10-26_v2\source_tables")
tab={t["sheet"]:t for t in json.loads((base/"supplementary_tables.json").read_text(encoding="utf-8"))["tables"]}
s8=tab["S8"]["rows"]; s9=tab["S9"]["rows"]; s12=tab["S12"]["rows"]
from collections import Counter
c8=Counter(r[1] for r in s8)
for r in s9:
    assert c8[r[0]]==r[2],(r[0],c8[r[0]],r[2])
assert sum(r[1] for r in s12)==9658 and sum(r[3] for r in s12)==9658
assert sum(r[9]=="PH" for r in s8)==1044
print("cross-table checks OK",c8,"S12 totals",sum(r[1] for r in s12),sum(r[3] for r in s12))
(Path(r"D:\桌面\single cell\修订交付_2026-10-26_v2\work")/"workbook_validation.txt").write_text("\n".join(issues) if issues else "PASS\n",encoding="utf-8")
