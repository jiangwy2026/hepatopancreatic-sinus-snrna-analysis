from pathlib import Path
import hashlib
root=Path(r"D:\桌面\single cell\github_release_2026-10-26")
skip={"MANIFEST_SHA256.txt","Figure5_source_panels","Figure6_source_panels"}
rows=[]
for p in sorted(root.rglob("*")):
    if not p.is_file(): continue
    if ".git" in p.parts or p.name in skip: continue
    rel=p.relative_to(root).as_posix()
    h=hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""): h.update(chunk)
    rows.append(f"{h.hexdigest()}  {rel}")
(root/"MANIFEST_SHA256.txt").write_text("\n".join(rows)+"\n",encoding="utf-8")
print("files",len(rows),"manifest",root/"MANIFEST_SHA256.txt")
