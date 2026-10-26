# -*- coding: utf-8 -*-
"""128_syncytium_verdict.py - 汇总：把被删的 PH 样群算进来后，PH 扩增倍数是否稳健"""
import io, sys, os
import numpy as np, pandas as pd
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
OUT = r"D:\桌面\single cell\analysis_output"
# 已知计数（来自 127 脚本）
kept = {"PH_blank": 366, "PH_treated": 678, "blank": 5349, "treated": 4309}
dropped = {"PH_blank": 126, "PH_treated": 206, "blank": 6396, "treated": 3732}
def frac(a): return 100 * a["PH_blank"] / a["blank"], 100 * a["PH_treated"] / a["treated"]
fb, ft = frac(kept); print("仅保留细胞：  blank %.2f%%  treated %.2f%%   -> %.2f 倍" % (fb, ft, ft / fb))
cb = {"PH_blank": kept["PH_blank"] + dropped["PH_blank"], "PH_treated": kept["PH_treated"] + dropped["PH_treated"],
      "blank": kept["blank"] + dropped["blank"], "treated": kept["treated"] + dropped["treated"]}
cb_b, cb_t = frac(cb); print("含被删 PH 样群：blank %.2f%%  treated %.2f%%   -> %.2f 倍" % (cb_b, cb_t, cb_t / cb_b))
print()
print("被删 PH 样群自身的扩增：blank %d 个 (%.2f%%) → treated %d 个 (%.2f%%)  -> %.2f 倍" % (
    dropped["PH_blank"], 100*dropped["PH_blank"]/dropped["blank"], dropped["PH_treated"],
    100*dropped["PH_treated"]/dropped["treated"], (dropped["PH_treated"]/dropped["treated"])/(dropped["PH_blank"]/dropped["blank"])))
