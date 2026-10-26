# -*- coding: utf-8 -*-
"""Generate an editable SVG model of the hepatopancreatic sinus hematopoietic niche."""

from pathlib import Path
from html import escape

ROOT = Path(r"D:\桌面\single cell")
OUT = ROOT / "analysis_output" / "figures_pub" / "2026-10-22_niche_mechanism"
ASSETS = OUT / "assets"
ASSETS.mkdir(parents=True, exist_ok=True)

W, H = 1720, 1190
parts = []
parts.append(f'''<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" role="img" aria-labelledby="title desc">
<title id="title">Hepatopancreatic sinus as an inducible hematopoietic niche</title>
<desc id="desc">Working model linking hemospasia-induced Hsp metabolic and protein remodeling, candidate ECM-integrin support from hepatic sinus niche cells, PH proliferation and expansion, and GH replenishment.</desc>
<defs>
  <style>
    text {{ font-family: 'Times New Roman', 'PingFang SC', 'Microsoft YaHei', 'Noto Sans CJK SC', system-ui, sans-serif; }}
    .title {{ font-size: 28px; font-weight: 700; fill: #0f172a; }}
    .subtitle {{ font-size: 16px; fill: #475569; }}
    .panel-title {{ font-size: 18px; font-weight: 700; fill: #334155; }}
    .box-title {{ font-size: 15px; font-weight: 700; fill: #1e293b; }}
    .body {{ font-size: 13px; fill: #334155; }}
    .small {{ font-size: 11px; fill: #64748b; }}
    .micro {{ font-size: 10px; fill: #64748b; }}
    .num {{ font-size: 22px; font-weight: 700; fill: #b2182b; }}
    .arrow-label {{ font-size: 11px; fill: #64748b; }}
    .group {{ fill: #f8fafc; stroke: #94a3b8; stroke-width: 1.2; stroke-dasharray: 7 5; }}
    .box {{ stroke-width: 1.5; }}
    .dash-sep {{ stroke: #cbd5e1; stroke-width: 1; stroke-dasharray: 4 4; }}
  </style>
  <marker id="arrow" markerWidth="8" markerHeight="8" refX="2" refY="4" orient="auto" markerUnits="userSpaceOnUse"><path d="M0,0 L8,4 L0,8 L2,4 z" fill="#64748b"/></marker>
  <marker id="arrow-blue" markerWidth="8" markerHeight="8" refX="2" refY="4" orient="auto" markerUnits="userSpaceOnUse"><path d="M0,0 L8,4 L0,8 L2,4 z" fill="#3b82f6"/></marker>
  <marker id="arrow-orange" markerWidth="8" markerHeight="8" refX="2" refY="4" orient="auto" markerUnits="userSpaceOnUse"><path d="M0,0 L8,4 L0,8 L2,4 z" fill="#f59e0b"/></marker>
  <marker id="arrow-green" markerWidth="8" markerHeight="8" refX="2" refY="4" orient="auto" markerUnits="userSpaceOnUse"><path d="M0,0 L8,4 L0,8 L2,4 z" fill="#22c55e"/></marker>
  <marker id="arrow-red" markerWidth="8" markerHeight="8" refX="2" refY="4" orient="auto" markerUnits="userSpaceOnUse"><path d="M0,0 L8,4 L0,8 L2,4 z" fill="#ef4444"/></marker>
  <marker id="arrow-purple" markerWidth="8" markerHeight="8" refX="2" refY="4" orient="auto" markerUnits="userSpaceOnUse"><path d="M0,0 L8,4 L0,8 L2,4 z" fill="#a855f7"/></marker>
</defs>
<rect x="0" y="0" width="1720" height="1190" fill="#ffffff"/>
''')


def text(x, y, value, cls="body", anchor="middle", weight=None, fill=None, extra=""):
    attrs = [f'x="{x}"', f'y="{y}"', f'class="{cls}"', f'text-anchor="{anchor}"']
    if weight is not None: attrs.append(f'font-weight="{weight}"')
    if fill is not None: attrs.append(f'fill="{fill}"')
    if extra: attrs.append(extra)
    return f'<text {" ".join(attrs)}>{escape(value)}</text>'


def lines(x, y, values, cls="body", gap=18, anchor="middle", weights=None, fills=None):
    out = []
    for i, value in enumerate(values):
        yy = y + i * gap
        weight = weights[i] if weights else None
        fill = fills[i] if fills else None
        out.append(text(x, yy, value, cls=cls, anchor=anchor, weight=weight, fill=fill))
    return "\n".join(out)


def rect(x, y, w, h, fill="#ffffff", stroke="#94a3b8", rx=10, cls="box", dash=None, sw=1.5):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}" stroke="{stroke}" class="{cls}" stroke-width="{sw}"{d}/>'


def path(d, color="#64748b", width=1.5, marker=None, dash=None, opacity=1):
    m = f' marker-end="url(#{marker})"' if marker else ""
    da = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<path d="{d}" fill="none" stroke="{color}" stroke-width="{width}" stroke-linecap="round" stroke-linejoin="round" opacity="{opacity}"{m}{da}/>'


def pill(x, y, w, label, fill, stroke, text_fill):
    return rect(x, y, w, 24, fill=fill, stroke=stroke, rx=12, cls="", sw=1) + text(x+w/2, y+16, label, cls="micro", fill=text_fill, weight="700")


# Title and framing
parts.append(text(860, 42, "Hepatopancreatic sinus as an inducible hematopoietic niche", cls="title"))
parts.append(text(860, 72, "Metabolic licensing and candidate ECM–integrin support converge on PH expansion after hemospasia", cls="subtitle"))
parts.append(rect(30, 105, 510, 755, cls="group", rx=12))
parts.append(rect(570, 105, 570, 755, cls="group", rx=12))
parts.append(rect(1170, 105, 520, 755, cls="group", rx=12))
parts.append(rect(30, 890, 1660, 255, cls="group", rx=12))

# Left: input and Hsp remodeling
parts.append(text(55, 135, "1. Hemospasia and Hsp remodeling", cls="panel-title", anchor="start"))
parts.append(rect(55, 160, 185, 70, fill="#fee2e2", stroke="#ef4444"))
parts.append(lines(147.5, 186, ["Hemospasia", "systemic hemocyte loss"], gap=20, weights=["700", None], fills=["#991b1b", "#7f1d1d"]))
parts.append(rect(285, 160, 225, 70, fill="#dbeafe", stroke="#3b82f6"))
parts.append(lines(397.5, 186, ["Hs plasma (Hsp)", "metabolic + protein cargo"], gap=20, weights=["700", None], fills=["#1e40af", "#1e3a8a"]))
parts.append(path("M245,195 L269,195", color="#3b82f6", marker="arrow-blue"))

# Hsp metabolic box
parts.append(path("M397.5,235 L397.5,269", color="#3b82f6", marker="arrow-blue"))
parts.append(rect(55, 285, 455, 125, fill="#eff6ff", stroke="#60a5fa"))
parts.append(text(282.5, 311, "Metabolic licensing", cls="box-title"))
parts.append(lines(282.5, 338, ["amino acids · nucleotide precursors · redox protection", "betaine supplementation increases THC", "growth-permissive Hsp environment"], gap=21, cls="body", fills=["#1e3a8a", "#1e40af", "#475569"]))

# Protein program
parts.append(path("M282.5,415 L282.5,439", color="#f59e0b", marker="arrow-orange"))
parts.append(rect(55, 450, 455, 180, fill="#fffbeb", stroke="#f59e0b"))
parts.append(text(282.5, 477, "Time-resolved protein program", cls="box-title"))
parts.append(lines(282.5, 504, ["6 h: eIF4B · ARHGAP12 · dynein light chain", "12–24 h: EF1α · RPS3A · RPS9", "sustained: filamin-A · β-parvin", "ECM/adhesion + biosynthetic maintenance"], gap=22, cls="body", fills=["#92400e", "#92400e", "#92400e", "#475569"]))
# syncytium box
parts.append(path("M282.5,635 L282.5,659", color="#f59e0b", marker="arrow-orange"))
parts.append(rect(55, 670, 455, 70, fill="#fff7ed", stroke="#fb923c"))
parts.append(lines(282.5, 696, ["Syncytium → mononuclear cells at 6 h", "morphological transition during recovery"], gap=20, weights=["700", None], fills=["#9a3412", "#7c2d12"]))
parts.append(text(282.5, 785, "Hsp data: 0 / 6 / 12 / 24 h proteome", cls="small"))
parts.append(text(282.5, 806, "Metabolome: 244 up / 39 down vs circulating plasma", cls="small"))

# Middle: niche cells
parts.append(text(595, 135, "2. Hs niche relay", cls="panel-title", anchor="start"))
parts.append(text(855, 165, "candidate niche cells", cls="small", weight="700"))
cell_specs = [
    (635, "MyoEpi", "#dcfce7", "#22c55e", "#166534"),
    (755, "M", "#f3e8ff", "#a855f7", "#6b21a8"),
    (875, "MH2", "#ffedd5", "#f97316", "#9a3412"),
    (995, "TS", "#dbeafe", "#3b82f6", "#1e40af"),
]
for cx, label, fill, stroke, txtfill in cell_specs:
    parts.append(f'<ellipse cx="{cx}" cy="195" rx="44" ry="29" fill="{fill}" stroke="{stroke}" stroke-width="1.4"/>')
    parts.append(f'<ellipse cx="{cx}" cy="195" rx="15" ry="11" fill="{stroke}" opacity="0.55"/>')
    parts.append(text(cx, 241, label, cls="body", weight="700", fill=txtfill))
parts.append(path("M875,252 L875,279", color="#22c55e", marker="arrow-green"))

# Ligand-receptor panel
parts.append(rect(595, 290, 520, 260, fill="#f0fdf4", stroke="#22c55e", rx=12))
parts.append(text(855, 316, "Candidate ECM / adhesion relay", cls="box-title"))
parts.append(pill(1002, 302, 96, "candidate", "#fef3c7", "#f59e0b", "#92400e"))
rows = [
    (335, "VWF/hemocytin → ITGB1", "MyoEpi → PH/GH   |   0.178 / 0.181", "#dcfce7", "#22c55e", "0.178 / 0.181"),
    (393, "THBS4 → ITGB1", "M/MH2 → PH/GH   |   0.048–0.070", "#ecfccb", "#84cc16", "0.048–0.070"),
    (451, "JAG1 → NOTCH1", "M → PH/GH   |   0.015–0.026", "#f1f5f9", "#94a3b8", "0.015–0.026"),
]
for y, pair, desc, fill, stroke, num in rows:
    parts.append(rect(615, y, 480, 48, fill=fill, stroke=stroke, rx=8, sw=1.2))
    parts.append(text(635, y+19, pair, cls="body", anchor="start", weight="700", fill="#1e293b"))
    parts.append(text(635, y+38, desc, cls="small", anchor="start", fill="#475569"))
parts.append(text(855, 526, "pathway sum: ECM–integrin 1.300 > Notch 0.063 > Wnt 0.00008", cls="small", fill="#166534"))
parts.append(text(855, 568, "receptor detection: PH 7.6% · GH 8.4% · niche 5.7%", cls="small"))
parts.append(text(855, 597, "low but detectable ITGB1; ligand and receptor evidence converge", cls="small"))
parts.append(text(855, 626, "No direct CXCL12/SCF orthologs: candidate relay restricted to", cls="small"))
parts.append(text(855, 646, "SonicParanoid-supported ECM–integrin / Notch components", cls="small"))
parts.append(pill(780, 680, 150, "functional validation pending", "#fef2f2", "#fca5a5", "#991b1b"))

# Curved arrow from relay to PH
parts.append(path("M1115,415 C1150,365 1280,270 1365,230", color="#22c55e", width=1.5, marker="arrow-green"))
parts.append(text(1265, 287, "candidate support", cls="small", fill="#166534"))

# Right: PH activation
parts.append(text(1195, 135, "3. PH activation and output", cls="panel-title", anchor="start"))
parts.append(f'<ellipse cx="1430" cy="205" rx="66" ry="48" fill="#fee2e2" stroke="#ef4444" stroke-width="2"/>')
parts.append(f'<ellipse cx="1430" cy="205" rx="24" ry="18" fill="#ef4444" opacity="0.72"/>')
parts.append(text(1430, 211, "PH", cls="box-title", weight="700", fill="#7f1d1d", extra='font-size="18"'))
parts.append(text(1430, 270, "prohemocyte", cls="small", fill="#991b1b", weight="700"))
parts.append(pill(1340, 282, 180, "BrdU+ increase in Hs", "#fee2e2", "#ef4444", "#991b1b"))

parts.append(rect(1200, 315, 460, 115, fill="#fff1f2", stroke="#fb7185", rx=10))
parts.append(text(1430, 341, "Proliferative response", cls="box-title"))
parts.append(lines(1430, 365, ["G2M 46.2% (n = 1,044)", "PH 6.8% → 15.7% (2.30×)", "dominant proliferative state"], gap=20, cls="body", fills=["#9f1239", "#9f1239", "#475569"]))

parts.append(rect(1200, 450, 460, 145, fill="#ffffff", stroke="#cbd5e1", rx=10))
parts.append(text(1430, 475, "PH substate remodeling", cls="box-title"))
substates = [
    (1220, 492, "s0 ECM/ribosome", "0.26×", "#e0f2fe", "#38bdf8", "#075985"),
    (1440, 492, "s1 progenitor", "1.54×", "#ffedd5", "#fb923c", "#9a3412"),
    (1220, 545, "s2 immune/ECM", "4.00×", "#fee2e2", "#f87171", "#991b1b"),
    (1440, 545, "s3 TGM1-low", "0.92×", "#dcfce7", "#4ade80", "#166534"),
]
for x, y, name, fc, fill, stroke, txtfill in substates:
    parts.append(rect(x, y, 200, 42, fill=fill, stroke=stroke, rx=8, sw=1.2))
    parts.append(text(x+12, y+18, name, cls="small", anchor="start", weight="700", fill=txtfill))
    parts.append(text(x+188, y+33, fc, cls="body", anchor="end", weight="700", fill=txtfill))

parts.append(rect(1200, 615, 460, 102, fill="#ecfdf5", stroke="#22c55e", rx=10))
parts.append(text(1430, 641, "Cycle-adjusted progenitor identity", cls="box-title", fill="#166534"))
parts.append(text(1265, 675, "z = 1.14", cls="num", anchor="start", fill="#15803d"))
parts.append(text(1265, 697, "91.1% of PH above expected", cls="body", anchor="start", fill="#166534"))
parts.append(text(1580, 675, "cycle-adjusted", cls="small", anchor="end", fill="#166534"))

# PH -> GH compact lineage
parts.append(f'<ellipse cx="1275" cy="775" rx="36" ry="26" fill="#fee2e2" stroke="#ef4444" stroke-width="1.5"/>')
parts.append(f'<ellipse cx="1275" cy="775" rx="12" ry="9" fill="#ef4444" opacity="0.7"/>')
parts.append(text(1275, 780, "PH", cls="micro", weight="700", fill="#7f1d1d"))
parts.append(path("M1315,775 L1368,775", color="#22c55e", marker="arrow-green"))
parts.append(text(1341, 759, "maturation", cls="micro", fill="#166534"))
parts.append(f'<ellipse cx="1465" cy="775" rx="40" ry="28" fill="#dcfce7" stroke="#22c55e" stroke-width="1.5"/>')
parts.append(f'<ellipse cx="1465" cy="775" rx="13" ry="9" fill="#22c55e" opacity="0.7"/>')
parts.append(text(1465, 780, "GH", cls="micro", weight="700", fill="#166534"))
parts.append(text(1430, 828, "PH → GH replenishment", cls="body", weight="700", fill="#166534"))
parts.append(text(1430, 848, "substate order: s1 → s0 → s2 → s3", cls="small"))

# Bottom working model
parts.append(text(55, 918, "Working model and cross-organ contrast", cls="panel-title", anchor="start"))
parts.append(rect(55, 945, 230, 80, fill="#eff6ff", stroke="#3b82f6", rx=10))
parts.append(lines(170, 973, ["Hsp metabolic licensing", "+ ECM–integrin support"], gap=20, weights=["700", None], fills=["#1e40af", "#1e3a8a"]))
parts.append(path("M290,985 L319,985", color="#64748b", marker="arrow"))
parts.append(rect(325, 945, 230, 80, fill="#fef3c7", stroke="#f59e0b", rx=10))
parts.append(lines(440, 973, ["PH cell-cycle entry", "+ PH expansion"], gap=20, weights=["700", None], fills=["#92400e", "#78350f"]))
parts.append(path("M560,985 L589,985", color="#64748b", marker="arrow"))
parts.append(rect(595, 945, 230, 80, fill="#dcfce7", stroke="#22c55e", rx=10))
parts.append(lines(710, 973, ["GH replenishment", "functional recovery"], gap=20, weights=["700", None], fills=["#166534", "#14532d"]))
parts.append(text(440, 1062, "in vivo BrdU response + THC recovery", cls="small", fill="#166534"))
parts.append(text(440, 1082, "metabolic and protein evidence are temporal, not causal", cls="micro"))

# Cross-organ contrast
parts.append(text(1110, 947, "Cross-organ distinction", cls="box-title", anchor="start", fill="#334155"))
parts.append(rect(1120, 970, 250, 100, fill="#fff1f2", stroke="#fb7185", rx=10))
parts.append(lines(1245, 995, ["Hs PH", "BrdU+ increase", "adjusted z = 1.14"], gap=20, weights=["700", None, "700"], fills=["#991b1b", "#9f1239", "#15803d"]))
parts.append(rect(1390, 970, 250, 100, fill="#eff6ff", stroke="#60a5fa", rx=10))
parts.append(lines(1515, 995, ["Circulating PH", "cycling program high", "BrdU+ not enriched"], gap=20, weights=["700", None, None], fills=["#1e40af", "#1d4ed8", "#475569"]))
parts.append(text(1380, 1100, "transcriptionally cycling ≠ injury-induced proliferation", cls="small", fill="#1e40af"))

# Footer note
parts.append(text(860, 1162, "Working model; candidate ECM/Notch axis awaits spatial and functional validation. Circulating PH is an external reference atlas, not a hemospasia time course.", cls="micro"))

parts.append("</svg>")
svg = "\n".join(parts)
out = ASSETS / "hs_niche_mechanism.svg"
out.write_text(svg, encoding="utf-8")
print(out)
