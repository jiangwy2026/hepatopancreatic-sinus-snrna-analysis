import pandas as pd, glob, os
OUT = r"D:\桌面\single cell\analysis_output"
for f in sorted(glob.glob(os.path.join(OUT, "elife_ref", "elife-80127-supp*-v2.csv"))):
    d = pd.read_csv(f, nrows=3)
    tc = d["Target_Cluster"].astype(str).tolist() if "Target_Cluster" in d.columns else "?"
    print(os.path.basename(f), "| Target_Cluster:", set(tc), "| cols:", list(d.columns)[:4])
