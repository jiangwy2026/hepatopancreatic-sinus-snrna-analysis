import os, sys, warnings
warnings.filterwarnings("ignore")
import pandas as pd
sys.stdout.reconfigure(encoding="utf-8")

ROOT = r"D:\桌面\single cell"
KEYS = ["human", "homo", "shrimp", "litopenaeus", "penaeus", "drosophila",
        "melanogaster", "danio", "rerio", "zebrafish", "ortholog", "orthogroup",
        "ortho", "fly", "ensp", "xp_", "species"]

def score_cols(cols):
    s = " | ".join(str(c).lower() for c in cols)
    return sum(1 for k in KEYS if k in s), s

hits = []
for dp, dn, fn in os.walk(ROOT):
    if any(x in dp for x in [".venv", "__pycache__", ".git", "env_r", "env_sonic", "mamba"]):
        continue
    for f in fn:
        lf = f.lower()
        p = os.path.join(dp, f)
        try:
            sz = os.path.getsize(p)
        except Exception:
            continue
        if sz < 500 or sz > 150_000_000:
            continue
        try:
            if lf.endswith((".xlsx", ".xls")):
                xl = pd.ExcelFile(p)
                for sh in xl.sheet_names:
                    d = pd.read_excel(p, sheet_name=sh, nrows=3)
                    sc, joined = score_cols(list(d.columns))
                    if sc >= 2:
                        hits.append((sc, sz, p, sh, len(d.columns), joined[:220]))
            elif lf.endswith((".csv", ".tsv", ".txt")):
                d = pd.read_csv(p, nrows=3, sep=None, engine="python", encoding="utf-8",
                                on_bad_lines="skip")
                sc, joined = score_cols(list(d.columns))
                if sc >= 2:
                    hits.append((sc, sz, p, "-", len(d.columns), joined[:220]))
        except Exception:
            continue

hits.sort(key=lambda x: (-x[0], -x[1]))
print(f"candidate tables: {len(hits)}")
for sc, sz, p, sh, nc, cols in hits:
    print(f"[{sc}] {sz:>12,d}  {p}  ::sheet={sh}  cols={nc}")
    print("      ", cols)
