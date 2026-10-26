import zipfile, re, os, io, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
root = r"D:\桌面\single cell\文章所需\cell report"
for fn in ("Manuscript_0702.docx", "MANUSCRIPTS_0702.docx", "Cover letter 0702.docx", "Cover letter 0703.docx"):
    p = os.path.join(root, fn)
    z = zipfile.ZipFile(p)
    names = z.namelist()
    xml = z.read("word/document.xml").decode("utf8", "replace")
    # paragraphs
    paras = re.findall(r"<w:p[ >].*?</w:p>", xml, flags=re.S)
    out = []
    for pa in paras:
        txt = "".join(re.findall(r"<w:t[^>]*>(.*?)</w:t>", pa, flags=re.S))
        txt = re.sub(r"<[^>]+>", "", txt)
        if txt.strip(): out.append(txt.strip())
    body = "\n".join(out)
    imgs = [n for n in names if n.startswith("word/media/")]
    print("=" * 90)
    print("%s | 段落 %d | 提取字符 %d | 内嵌图片 %d (%s)" % (fn, len(out), len(body), len(imgs),
          ", ".join(sorted({os.path.splitext(i)[1] for i in imgs}))))
    print("=" * 90)
    print(body[:3000])
    print("...\n")
    open(os.path.join(os.path.dirname(p), "_extract_" + fn.replace(".docx", ".txt")), "w", encoding="utf8").write(body)
