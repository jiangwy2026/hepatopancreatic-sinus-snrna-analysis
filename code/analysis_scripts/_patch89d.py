p = r"D:\桌面\single cell\analysis_output\scripts\89_velocity_scvelo.py"
src = open(p, encoding="utf8").read()
old = "import os, sys, json\nimport numpy as np, pandas as pd"
new = "import os, sys, json\nimport matplotlib\nmatplotlib.use(\"Agg\")          # must be set before pyplot/scvelo import (no Tk in this environment)\nimport numpy as np, pandas as pd"
assert old in src
src = src.replace(old, new, 1)
src = src.replace("    import matplotlib; matplotlib.use(\"Agg\")\n", "")
open(p, "w", encoding="utf8").write(src)
print("patched")
