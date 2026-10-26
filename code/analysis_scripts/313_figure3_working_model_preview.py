# -*- coding: utf-8 -*-
"""Generate a preview of the Figure 3 working model with evidence-strength arrows."""

from pathlib import Path
from html import escape

ROOT = Path(r"D:\桌面\single cell")
OUT = ROOT / "analysis_output" / "figures_pub" / "2026-10-23_Figure3_working_model_preview"
OUT.mkdir(parents=True, exist_ok=True)
W, H = 1600, 920
parts = []
parts.append(f'''<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">
<title>Working model for metabolic–ECM licensing of Hs PH expansion</title>
<defs>
<style>text {{ font-family: 'Times New Roman', 'PingFang SC', 'Microsoft YaHei', 'Noto Sans CJK SC', sans-serif; }} .title {{ font-size: 27px; font-weight: 700; fill:#0f172a; }} .subtitle {{ font-size: 13px; fill:#64748b; }} .panel {{ font-size: 17px; font-weight: 700; fill:#334155; }} .box {{ font-size: 14px; font-weight: 700; fill:#1e293b; }} .body {{ font-size: 12px; fill:#334155; }} .small {{ font-size: 10px; fill:#64748b; }} .big {{ font-size: 23px; font-weight: 700; fill:#b2182b; }} .group {{ fill:#f8fafc; stroke:#94a3b8; stroke-width:1.2; stroke-dasharray:7 5; }} </style>
<marker id="arrow" markerWidth="8" markerHeight="8" refX="2" refY="4" orient="auto" markerUnits="userSpaceOnUse"><path d="M0,0 L8,4 L0,8 L2,4 z" fill="#64748b"/></marker>
<marker id="arrow-red" markerWidth="8" markerHeight="8" refX="2" refY="4" orient="auto" markerUnits="userSpaceOnUse"><path d="M0,0 L8,4 L0,8 L2,4 z" fill="#ef4444"/></marker>
<marker id="arrow-green" markerWidth="8" markerHeight="8" refX="2" refY="4" orient="auto" markerUnits="userSpaceOnUse"><path d="M0,0 L8,4 L0,8 L2,4 z" fill="#22c55e"/></marker>
</defs>
<rect x="0" y="0" width="1600" height="920" fill="white"/>
''')


def text(x, y, value, cls='body', anchor='middle', weight=None, fill=None, size=None):
    attrs=[f'x="{x}"',f'y="{y}"',f'class="{cls}"',f'text-anchor="{anchor}"']
    if weight: attrs.append(f'font-weight="{weight}"')
    if fill: attrs.append(f'fill="{fill}"')
    if size: attrs.append(f'font-size="{size}"')
    return f'<text {" ".join(attrs)}>{escape(value)}</text>'


def lines(x,y,vals,gap=18,cls='body',weights=None,fills=None):
    return '\n'.join(text(x,y+i*gap,v,cls=cls,weight=(weights[i] if weights else None),fill=(fills[i] if fills else None)) for i,v in enumerate(vals))


def rect(x,y,w,h,fill='white',stroke='#94a3b8',rx=10,sw=1.4,dash=None,cls=''):
    d=f' stroke-dasharray="{dash}"' if dash else ''
    c=f' class="{cls}"' if cls else ''
    return f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"{c}{d}/>'


def path(d,color='#64748b',width=1.5,marker='arrow',dash=None):
    da=f' stroke-dasharray="{dash}"' if dash else ''
    return f'<path d="{d}" fill="none" stroke="{color}" stroke-width="{width}" stroke-linecap="round" stroke-linejoin="round" marker-end="url(#{marker})"{da}/>'


parts.append(text(800,42,'A working model for metabolic–ECM licensing of Hs PH expansion','title'))
parts.append(text(800,70,'Solid arrows: observed association    Dashed arrows: candidate or untested link','subtitle'))
parts.append(rect(30,105,430,590,cls='group',rx=12))
parts.append(rect(485,105,500,590,cls='group',rx=12))
parts.append(rect(1010,105,560,590,cls='group',rx=12))
parts.append(rect(30,725,1540,165,cls='group',rx=12))

# Left: observed input
parts.append(text(55,135,'1. Hemospasia and Hsp remodeling','panel',anchor='start'))
parts.append(rect(55,165,150,58,fill='#fee2e2',stroke='#ef4444'))
parts.append(lines(130,194,['Hemospasia','hemocyte loss'],gap=18,weights=['700',None],fills=['#991b1b','#7f1d1d']))
parts.append(rect(250,165,180,58,fill='#dbeafe',stroke='#3b82f6'))
parts.append(lines(340,194,['Hs plasma (Hsp)','metabolic + protein'],gap=18,weights=['700',None],fills=['#1e40af','#1e3a8a']))
parts.append(path('M205,194 L239,194',color='#3b82f6'))
parts.append(path('M340,223 L340,258',color='#3b82f6'))
parts.append(rect(55,270,375,130,fill='#eff6ff',stroke='#60a5fa'))
parts.append(text(242,297,'Metabolic licensing','box'))
parts.append(lines(242,323,['amino acids · nucleotides · redox','betaine supports THC recovery','growth-permissive Hsp'],gap=21,fills=['#1e3a8a','#1e40af','#475569']))
parts.append(path('M242,400 L242,435',color='#f59e0b'))
parts.append(rect(55,445,375,145,fill='#fffbeb',stroke='#f59e0b'))
parts.append(text(242,472,'Protein program','box'))
parts.append(lines(242,498,['6 h: eIF4B · ARHGAP12 · dynein light chain','12–24 h: EF1α · RPS3A · RPS9','sustained: filamin-A · β-parvin'],gap=21,fills=['#92400e','#92400e','#475569']))
parts.append(text(242,630,'Measured in Hsp','small',weight='700'))
parts.append(text(242,650,'metabolome and time-resolved proteome','small'))

# Middle: candidate niche
parts.append(text(735,135,'2. Candidate Hs niche','panel',anchor='start'))
parts.append(text(735,165,'candidate niche cells','small',weight='700'))
for cx,label,fill,stroke,filltext in [(600,'MyoEpi','#dcfce7','#22c55e','#166534'),(735,'M','#f3e8ff','#a855f7','#6b21a8'),(870,'MH2','#ffedd5','#f97316','#9a3412')]:
    parts.append(f'<ellipse cx="{cx}" cy="200" rx="48" ry="31" fill="{fill}" stroke="{stroke}" stroke-width="1.5"/>')
    parts.append(f'<ellipse cx="{cx}" cy="200" rx="16" ry="12" fill="{stroke}" opacity="0.6"/>')
    parts.append(text(cx,245,label,'body',weight='700',fill=filltext))
parts.append(path('M735,260 L735,300',color='#22c55e',marker='arrow-green',dash='8 6'))
parts.append(rect(535,305,400,220,fill='#f0fdf4',stroke='#22c55e'))
parts.append(text(735,332,'Candidate ECM / adhesion relay','box'))
parts.append(rect(555,350,360,46,fill='#dcfce7',stroke='#22c55e',rx=7,sw=1))
parts.append(text(735,370,'VWF/hemocytin → ITGB1','body',weight='700'))
parts.append(text(735,387,'MyoEpi → PH/GH','small'))
parts.append(rect(555,407,360,46,fill='#ecfccb',stroke='#84cc16',rx=7,sw=1))
parts.append(text(735,427,'THBS4 → ITGB1','body',weight='700'))
parts.append(text(735,444,'M/MH2 → PH/GH','small'))
parts.append(rect(555,464,360,40,fill='#f1f5f9',stroke='#94a3b8',rx=7,sw=1))
parts.append(text(735,488,'JAG1 → NOTCH1  (weak)','small'))
parts.append(text(735,548,'ITGB1 detection: PH 7.6% · GH 8.4%','small'))
parts.append(text(735,568,'candidate axis awaits spatial/functional validation','small',fill='#64748b'))
# dashed arrow to PH
parts.append(path('M935,415 C965,360 1040,300 1080,250',color='#22c55e',marker='arrow-green',dash='8 6'))
parts.append(text(1010,320,'candidate','small',fill='#166534'))

# Right: observed PH response
parts.append(text(1040,135,'3. PH response','panel',anchor='start'))
parts.append(f'<ellipse cx="1290" cy="205" rx="70" ry="50" fill="#fee2e2" stroke="#ef4444" stroke-width="2"/>')
parts.append(f'<ellipse cx="1290" cy="205" rx="24" ry="18" fill="#ef4444" opacity="0.7"/>')
parts.append(text(1290,211,'PH','box',weight='700',fill='#7f1d1d',size='18'))
parts.append(text(1290,268,'prohemocyte','small',weight='700',fill='#991b1b'))
parts.append(rect(1040,300,480,92,fill='#fff1f2',stroke='#fb7185'))
parts.append(text(1280,326,'Observed PH activation','box'))
parts.append(lines(1280,350,['G2M 46.2%  ·  PH 6.8% → 15.7% (2.30×)','s1/s2 expansion  ·  BrdU⁺ PH increases in Hs'],gap=20,fills=['#9f1239','#9f1239']))
parts.append(rect(1040,415,480,95,fill='#ecfdf5',stroke='#22c55e'))
parts.append(text(1280,441,'Cycle-adjusted progenitor identity','box',fill='#166534'))
parts.append(text(1090,478,'z = 1.14','big',anchor='start',fill='#15803d'))
parts.append(text(1300,478,'91.1% above expected','body',anchor='start',fill='#166534'))
# PH -> GH dashed
parts.append(f'<ellipse cx="1110" cy="590" rx="38" ry="27" fill="#fee2e2" stroke="#ef4444"/>')
parts.append(text(1110,595,'PH','small',weight='700',fill='#7f1d1d'))
parts.append(path('M1150,590 L1198,590',color='#64748b',marker='arrow',dash='8 6'))
parts.append(text(1174,570,'maturation','small',fill='#64748b'))
parts.append(f'<ellipse cx="1245" cy="590" rx="40" ry="28" fill="#dcfce7" stroke="#22c55e"/>')
parts.append(text(1245,595,'GH','small',weight='700',fill='#166534'))
parts.append(text(1280,626,'PH → GH: trajectory-supported, not lineage-traced','small',fill='#64748b'))

# Bottom: evidence boundary
parts.append(text(55,752,'Evidence and claim boundary','box',anchor='start',fill='#334155'))
parts.append(rect(55,775,650,85,fill='#ffffff',stroke='#cbd5e1',rx=9,sw=1))
parts.append(text(75,798,'Observed:', 'body', anchor='start', weight='700', fill='#166534'))
parts.append(text(150,798,'BrdU⁺ Hs increase · Hsp metabolome/proteome · PH expansion','body',anchor='start'))
parts.append(text(75,826,'Candidate:', 'body', anchor='start', weight='700', fill='#92400e'))
parts.append(text(150,826,'ECM–integrin/Notch relay · niche-cell support · PH→GH direction','body',anchor='start'))
parts.append(rect(735,775,810,85,fill='#fff7ed',stroke='#fb923c',rx=9,sw=1))
parts.append(text(760,798,'Working model:', 'body', anchor='start', weight='700', fill='#9a3412'))
parts.append(text(890,798,'Hsp metabolic licensing + ECM–integrin support may create a proliferation-permissive Hs niche.','body',anchor='start'))
parts.append(text(760,826,'Direct niche function requires spatial co-localization and perturbation experiments.','body',anchor='start',fill='#9a3412'))
parts.append('</svg>')
out=OUT/'Figure3_working_model_preview.svg'
out.write_text('\n'.join(parts),encoding='utf-8')
print(out)
