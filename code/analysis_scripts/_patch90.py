p = r"D:\桌面\single cell\analysis_output\scripts\90_scenic_finish.py"
src = open(p, encoding="utf8").read()
old = "import numpy as np, pandas as pd\n"
new = """import numpy as np, pandas as pd
if not hasattr(np, "object"):      # numpy>=1.24 removed the alias pyscenic still uses
    np.object = object
import dask
dask.config.set({"dataframe.query-planning": False})   # legacy dask.dataframe API without dask-expr
"""
assert old in src
src = src.replace(old, new, 1)
src = src.replace('"""90_scenic_finish.py - SCENIC cisTarget pruning + AUCell (E:\\workspace\\scenic_env\\python.exe)."""',
                  'r"""90_scenic_finish.py - SCENIC cisTarget pruning + AUCell (E:\\workspace\\scenic_env\\python.exe)."""')
open(p, "w", encoding="utf8").write(src)
print("patched")
