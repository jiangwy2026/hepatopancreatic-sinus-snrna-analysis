# -*- coding: utf-8 -*-
"""191_individual_violin_exports.py
Render each volcano-labelled gene as a standalone violin and export SVG + 1200-dpi TIFF.
"""
import os, sys, re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
from scipy.stats import mannwhitneyu
from PIL import Image
sys.path.insert(0, r"D:\桌面\single cell\analysis_output")
from pub_style import apply_pub_style
apply_pub_style(7)
sc.settings.verbosity = 0
ROOT = r"D:\桌面\single cell"
OUT = os.path.join(ROOT, "analysis_output")
PREV = os.path.join(OUT, "figures_pub", "sc_summary_new_names", "individual_previews", "violin_single")
DELIV = os.path.join(ROOT, "文章配图", "最终_2026-10-17", "机制_单基因小提琴")
os.makedirs(PREV, exist_ok=True); os.makedirs(DELIV, exist_ok=True)
CTRL = "#4C78A8"; TRT = "#D6604D"; DPI = 1200
GENES = [
    ("ncbi_113823945", "hemocyte TGM", "pos"),
    ("ncbi_113814652", "ccnd2", "pos"),
    ("ncbi_113824171", "SSPO", "pos"),
    ("ncbi_113808957", "Notch-like", "pos"),
    ("ncbi_113824160", "THBS4", "neg"),
    ("ncbi_113808997", "PEN-3", "neg"),
    ("ncbi_113801825", "crustin", "neg"),
    ("ncbi_113810905", "WAP", "neg"),
]

def save_pub(fig, stem):
    svg = stem + ".svg"; tif = stem + ".tiff"
    fig.savefig(svg, bbox_inches="tight", pad_inches=0.04, facecolor="white")
    fig.savefig(tif, dpi=DPI, bbox_inches="tight", pad_inches=0.04, facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    with Image.open(tif) as im:
        if "A" in im.getbands():
            rgba = im.convert("RGBA"); bg = Image.new("RGB", rgba.size, (255,255,255))
            bg.paste(rgba, mask=rgba.getchannel("A")); tmp = tif + ".tmp.tiff"
            bg.save(tmp, format="TIFF", dpi=(DPI,DPI), compression="tiff_lzw")
            os.replace(tmp, tif)
    print("saved", svg); print("saved", tif)


def main():
    A = ad.read_h5ad(os.path.join(OUT, "adata_C2_within.h5ad"))
    A.X = A.layers["counts"].copy(); sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
    A.obs["pt_bin3"] = pd.qcut(A.obs["dpt_pseudotime"].astype(float), 3, labels=["P1","P2","P3"])
    E = A[A.obs["pt_bin3"] == "P1"].copy()
    gid2i = {str(g): i for i,g in enumerate(A.var_names.astype(str))}
    grp = E.obs["group"].astype(str).values
    rng = np.random.default_rng(3)
    rows = []
    for idx, (gid, label, kind) in enumerate(GENES, start=1):
        v = E[:, gid2i[gid]].X.toarray().ravel().astype(float)
        data = [v[grp == "blank"], v[grp == "treated"]]
        fig = plt.figure(figsize=(2.9, 3.35), facecolor="white")
        ax = fig.add_axes([0.20, 0.20, 0.76, 0.68])
        vp = ax.violinplot(data, positions=[1,2], widths=0.62, showextrema=False, showmeans=False)
        for body,col in zip(vp["bodies"], (CTRL,TRT)):
            body.set_facecolor(col); body.set_edgecolor("none"); body.set_alpha(0.55)
        for pos,val,col in zip((1,2), data, (CTRL,TRT)):
            j = rng.uniform(-0.07,0.07,size=len(val))
            ax.scatter(pos+j, val, s=0.9, c=col, alpha=0.42, linewidths=0, rasterized=True)
            ax.scatter([pos],[np.median(val)],s=11,c="white",edgecolors=col,linewidths=0.6,zorder=5)
        p = mannwhitneyu(data[1],data[0],alternative="two-sided").pvalue
        ymax=max(np.max(data[0]),np.max(data[1])); top=ymax+(0.25 if ymax>0 else 0.1)
        ax.set_ylim(-0.08, top*1.08+0.05)
        ax.plot([1,2],[top,top],color="#333333",lw=0.45)
        ax.text(1.5, top+0.025*top, "P = %.1e" % p, ha="center", va="bottom", fontsize=5.0)
        ax.set_xticks([1,2]); ax.set_xticklabels(["Control","Hemospasia"], fontsize=5.0, rotation=20, ha="right")
        ax.set_ylabel("Log-normalized expression", fontsize=5.3)
        ax.tick_params(axis="y",labelsize=4.8,length=1.1,width=0.45); ax.tick_params(axis="x",length=0,pad=1)
        for sp in ("top","right"): ax.spines[sp].set_visible(False)
        col = "#B2182B" if kind=="pos" else "#2166AC"
        ax.set_title(label, fontsize=6.4, color=col, pad=5)
        safe = re.sub(r"[^A-Za-z0-9._-]+","_",label)
        stem = os.path.join(DELIV, "%02d_%s" % (idx, safe))
        save_pub(fig, stem)
        png = os.path.join(PREV, "%02d_%s_preview.png" % (idx, safe))
        fig.savefig(png, dpi=300, bbox_inches="tight", pad_inches=0.04, facecolor="white")
        plt.close(fig)
        rows.append({"index":idx,"gene":label,"gid":gid,"class":kind,"p":p,
                     "svg":os.path.basename(stem)+".svg","tiff":os.path.basename(stem)+".tiff"})
    pd.DataFrame(rows).to_csv(os.path.join(OUT,"机制_单基因小提琴_manifest.csv"),index=False,encoding="utf-8-sig")
    print("DONE")

if __name__ == "__main__":
    main()
