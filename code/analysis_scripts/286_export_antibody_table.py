# -*- coding: utf-8 -*-
import os
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
from openpyxl.utils import get_column_letter

out = r"D:\桌面\肝血窦PH_marker抗体询价表_2026-10-21.xlsx"
rows = [
 ["1 首选","GRIK2","113823802","Atlas Antibodies","HPA014395","兔多抗","GRIK2","44.7%","89%","aa 20–880","62.3%","11.4%","首选先试","询价时要求提供 immunogen 序列，确认落在 aa 20–880"],
 ["2 备选","GRIK2","113823802","Atlas Antibodies","HPA014623","兔多抗","GRIK2","44.7%","89%","aa 20–880","62.3%","11.4%","备选","同上"],
 ["3 第二选择","KCNMA1 (SLO)","113807409","Atlas Antibodies","HPA054648","兔多抗","KCNMA1","59.8%","98%","aa 226–1154","65.9%","21.2%","第二选择","保守性最好，但其他群也表达；需与 GRIK2 共染"],
 ["4 增殖共标","CCND3","113814652","Atlas Antibodies","HPA042871","兔多抗","CCND3","41.9%","93%","aa 2–265","50.8%","9.2%","增殖共标","注意：对虾 PH 表达的是 CCND3-like，不是 CCND2"],
 ["5 备选","CCND3","113814652","Atlas Antibodies","CAB000116","兔多抗","CCND3","41.9%","93%","aa 2–265","50.8%","9.2%","备选","同上"],
 ["6 备选","TGM1","113823934","Atlas Antibodies","HPA015278","兔多抗","TGM1","38.1%","67%","aa 125–789","76.2%","49.6%","备选","PH 特异性一般；先确认 immunogen 在 aa 125–789"],
 ["7 建议定制","TGM4-like","113823945","—","定制","兔多抗","TGM4","39.0%","59%","aa 45–119","58.7%","10.2%","建议定制虾肽抗体","保守段太小；商业 TGM4 抗体大概率不交叉"],
 ["8 谨慎","TMPRSS9 (PPAF2)","113812814","Atlas Antibodies","HPA051483","兔多抗","TMPRSS9","34.4%","50%","aa 220–467","76.0%","21.5%","谨慎","保守段小，交叉反应不确定"],
 ["9 备选","SERPINB2 (SERPINB1)","113799939","Atlas Antibodies","HPA015480","兔多抗","SERPINB2","32.8%","83%","aa 10–415","62.8%","27.4%","备选","全长保守，但 PH 特异性一般"],
]
header = ["优先级","靶点","虾 GeneID","公司","货号","抗体类型","人同源物","Identity","Coverage","人保守区 (aa)","PH 阳性率","其他群阳性率","建议","备注"]
wb = Workbook(); ws = wb.active; ws.title = "抗体询价表"
ws.append(header)
for r in rows: ws.append(r)
thin = Side(style="thin", color="D9D9D9")
fill = PatternFill("solid", fgColor="1F4E79")
for c in range(1, len(header)+1):
    cell = ws.cell(row=1, column=c)
    cell.font = Font(name="Times New Roman", bold=True, color="FFFFFF", size=10)
    cell.fill = fill
    cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
    cell.border = Border(bottom=thin)
for row in ws.iter_rows(min_row=2, max_row=ws.max_row):
    for cell in row:
        cell.font = Font(name="Times New Roman", size=9.5)
        cell.alignment = Alignment(vertical="center", wrap_text=True)
        cell.border = Border(bottom=thin)
widths = [12,18,12,18,14,10,14,10,10,14,11,12,22,46]
for i,w in enumerate(widths,1): ws.column_dimensions[get_column_letter(i)].width = w
ws.freeze_panes = "A2"
ws.auto_filter.ref = f"A1:N{ws.max_row}"
# note
ws.append([])
n = ws.max_row + 1
ws.cell(row=n, column=1, value="说明：HPA/Atlas 抗体货号来自 Human Protein Atlas；HPA 的公开 API 不直接给出 immunogen 氨基酸范围，下单前必须让 Atlas 提供 immunogen 序列并核对是否落在人保守区内。CST #3741 为 anti-Cyclin D2，不推荐用于本项目的 PH 标记（对虾 PH 表达的是 CCND3-like）。")
ws.merge_cells(start_row=n, start_column=1, end_row=n, end_column=14)
ws.cell(row=n, column=1).font = Font(name="Times New Roman", size=9, italic=True)
ws.cell(row=n, column=1).alignment = Alignment(wrap_text=True, vertical="top")
ws.row_dimensions[n].height = 55
wb.save(out)
print("saved", out, os.path.getsize(out))
