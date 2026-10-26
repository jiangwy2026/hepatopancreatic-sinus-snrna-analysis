import pandas as pd, os, sys
sys.stdout.reconfigure(encoding="utf-8")
files = [
 r"D:\桌面\single cell\文章所需\Supplemental Table 5.xlsx",
 r"D:\桌面\single cell\文章所需\Supplemental Table 6.xlsx",
 r"D:\桌面\single cell\文章所需\PANS\肝血窦数据\Supplemental Table 6.xlsx",
 r"D:\桌面\single cell\文章所需\PANS\肝血窦数据\Supplemental Table 7.xlsx",
 r"D:\桌面\single cell\文章所需\图5数据\肝胰腺比对材料.xlsx",
 r"D:\桌面\single cell\文章所需\图5数据\淋巴比对材料.xlsx",
]
for p in files:
    print("="*70)
    print(p, os.path.getsize(p))
    try:
        xl = pd.ExcelFile(p)
        print("sheets:", xl.sheet_names)
        for sh in xl.sheet_names:
            d = pd.read_excel(p, sheet_name=sh)
            print(f"  {sh}: shape={d.shape} cols={list(d.columns)[:6]}")
    except Exception as e:
        print("ERR", e)
