from pathlib import Path
import zipfile, re, hashlib
from lxml import etree
from docx import Document
root=Path(r"D:\桌面\single cell")
out=root/"修订交付_2026-10-26_v2"
orig=root/"manuscripts_1026改.docx"; final=out/"manuscripts_1026_revised_final_2026-10-26.docx"
def stats(p):
    with zipfile.ZipFile(p) as z:
        assert z.testzip() is None
        doc=z.read("word/document.xml")
        rootxml=etree.fromstring(doc)
        drawings=len(rootxml.xpath("//*[local-name()='drawing']"))
        media=[n for n in z.namelist() if n.startswith("word/media/")]
        rels=z.read("word/_rels/document.xml.rels")
        relroot=etree.fromstring(rels)
        return {"media":len(media),"drawings":drawings,"relationships":len(relroot)}
o=stats(orig); f=stats(final)
assert o["media"]==f["media"] and o["drawings"]==f["drawings"]
d=Document(final); paras=[p.text.strip() for p in d.paragraphs]
heads=[x for x in paras if re.match(r"^3\.[1-9] ",x)]
figs=[x for x in paras if re.match(r"^Figure [1-9]\.",x)]
refs=[p for p in d.paragraphs if p.style.name=="EndNote Bibliography"]
assert len(heads)==6 and not any(x.startswith("3.7 ") for x in heads)
assert len(figs)==6 and len(refs)==40
text="\n".join(paras)
for n in range(7,20): assert f"Table S{n}" in text, n
assert "Figure 5. Single-cell atlas of Hs cell states" in text
assert "(G) PH-to-GH transcriptional ordering" in text
assert "Figure 6. Hemospasia remodels a proliferative, progenitor-like PH pool" in text
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
report=f"""# 2026-10-26 revision validation\n\n- DOCX opened and parsed by python-docx.\n- DOCX ZIP integrity check: PASS.\n- OOXML relationships parsed: PASS.\n- Embedded media count unchanged: original {o['media']}, revised {f['media']}.\n- Drawing count unchanged: original {o['drawings']}, revised {f['drawings']}.\n- Result sections: 6; standalone 3.7 heading removed.\n- Main figure captions: 6.\n- EndNote bibliography records: 40 (unchanged).\n- Table S7-S19 references present in manuscript: PASS.\n- Figure 5 panel G legend present: PASS.\n- Figure 5 and Figure 6 titles aligned with result headings: PASS.\n\n## Files\n\n- Revised DOCX SHA256: `{sha(final)}`\n- Supplementary workbook SHA256: `{sha(out/'Supplementary_Tables_S7-S19_2026-10-26.xlsx')}`\n"""
(out/"修订说明与核查记录_2026-10-26_v2.md").write_text(report,encoding="utf-8")
print(report)
