import re, sys, os, subprocess
sys.stdout.reconfigure(encoding="utf-8")
FA = r"E:\workspace\virtualbox_share\input_proteomes\pvannamei.faa"
want = ["XP_069984471","XP_069989473","XP_069992807","XP_069977048","XP_069980464",
        "XP_069985836","XP_069991869","XP_069971667","XP_069987369","XP_069969759",
        "XP_027231548","XP_069989416","XP_069989551","XP_027225409","XP_069998165",
        "XP_069988451","XP_069969627","XP_069986831","XP_027239322","XP_069973256",
        "XP_069973257","XP_069983458","XP_069974570","XP_069986011"]
out = []
keep = False
hdr = ""
for line in open(FA, encoding="utf-8", errors="ignore"):
    if line.startswith(">"):
        acc = line[1:].split()[0].split(".")[0]
        keep = acc in want
        hdr = line
        if keep:
            print("HEADER:", line.strip()[:140])
            out.append(line)
    elif keep:
        out.append(line)
print("found", sum(1 for l in out if l.startswith(">")), "of", len(want))
p = os.path.join(os.environ["TEMP"], "blastwork", "shrimp_new_lr.faa")
open(p, "w", encoding="utf-8").write("".join(out))
print(p)
