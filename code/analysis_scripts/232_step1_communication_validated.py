# -*- coding: utf-8 -*-
"""231_step1_communication.py

Step 1b: CellChat-style ligand-receptor communication from candidate niche cells
to hemocyte targets (PH, GH).

R/CellChat is not installed in this environment, and CellChatDB is a
human/mouse database that is not valid for Litopenaeus vannamei. This script
therefore implements the CellChat probability model (mean ligand x mean
receptor) with a label-permutation null on a curated shrimp LR panel built from
NCBI product-name assignments.
"""
import os, json, itertools
import numpy as np, pandas as pd, anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, FancyArrowPatch
from matplotlib.path import Path
from matplotlib.patches import PathPatch

OUT = r"D:\桌面\single cell\analysis_output"
STEP = os.path.join(OUT, "step1")
os.makedirs(STEP, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", os.path.join(OUT, ".mplconfig"))
RNG = np.random.default_rng(42)

SOURCES = ["TS", "MH2", "MyoEpi", "M"]
TARGETS = ["PH", "GH"]
ALL = SOURCES + TARGETS

# ligand, receptor, axis, product note
LR = [
 ("113821858","113802681","Notch","JAG1 -> NOTCH1"),
 ("113821858","113802400","Notch","JAG1 -> NOTCH2"),
 ("113815864","113802681","Notch","DLL1 -> NOTCH1"),
 ("113815864","113802400","Notch","DLL1 -> NOTCH2"),
 ("113806005","113826277","Wnt","WNT5A -> FZD5"),
 ("113806005","113826275","Wnt","WNT5A -> FZD4"),
 ("113804527","113823485","Angpt-Tie","ANGPT2 -> TEK-like"),
 ("113823118","113821319","ECM-integrin","THBS4 -> ITGAV"),
 ("113823118","113830322","ECM-integrin","THBS4 -> ITGB1"),
 ("113802906","113821319","ECM-integrin","LAMA1 -> ITGAV"),
 ("113802906","113830322","ECM-integrin","LAMA1 -> ITGB1"),
 ("113802906","113827435","ECM-integrin","LAMA1 -> ITGA4"),
 ("113805211","113821319","ECM-integrin","POSTN -> ITGAV"),
 ("113824171","113821319","ECM-integrin","VWF/hemocytin -> ITGAV"),
]


def load():
    a = ad.read_h5ad(os.path.join(STEP, "step1_niche.h5ad"))
    a.X = a.layers["lognorm"].copy()
    gid2var = {str(g).replace("ncbi_", ""): v for g, v in zip(a.var["gid"], a.var_names)}
    return a, gid2var


def means(a, idx, groups):
    X = np.asarray(a[:, idx].X.todense()).ravel()
    out = {}
    for g in groups:
        m = (a.obs["cell_type"].astype(str) == g).values
        out[g] = (float(X[m].mean()) if m.sum() else 0.0, float((X[m] > 0).mean()) if m.sum() else 0.0)
    return out


def main():
    a, g2v = load()
    rows = []
    n_perm = 1000
    # precompute per-pair cell vectors
    cache = {}
    for L, R, axis, note in LR:
        if L not in g2v or R not in g2v:
            continue
        cache[(L, R)] = (np.asarray(a[:, g2v[L]].X.todense()).ravel(),
                         np.asarray(a[:, g2v[R]].X.todense()).ravel())
    labels = a.obs["cell_type"].astype(str).values
    keep = np.isin(labels, ALL)
    lab = labels[keep]

    for (L, R), (xl, xr) in cache.items():
        axis = [x[2] for x in LR if x[0] == L and x[1] == R][0]
        note = [x[3] for x in LR if x[0] == L and x[1] == R][0]
        ml = means(a, g2v[L], ALL); mr = means(a, g2v[R], ALL)
        xl_k, xr_k = xl[keep], xr[keep]
        for s in SOURCES:
            for t in TARGETS:
                if s == t: continue
                ms, ts = lab == s, lab == t
                prob = ml[s][0] * mr[t][0]
                # permutation null: shuffle labels of the two cell sets
                null = np.empty(n_perm)
                for k in range(n_perm):
                    perm = RNG.permutation(lab)
                    null[k] = xl_k[perm == s].mean() * xr_k[perm == t].mean()
                p = (np.sum(null >= prob) + 1) / (n_perm + 1)
                rows.append({"ligand_gid":L,"receptor_gid":R,"axis":axis,"pair":note,
                             "source":s,"target":t,
                             "ligand_mean":ml[s][0],"ligand_pct":ml[s][1],
                             "receptor_mean":mr[t][0],"receptor_pct":mr[t][1],
                             "prob":prob,"p_perm":p,"n_source":int(ms.sum()),"n_target":int(ts.sum())})
    df = pd.DataFrame(rows).sort_values(["prob"], ascending=False)
    df.to_csv(os.path.join(STEP, "step1_lr_communication_validated.csv"), index=False, encoding="utf-8-sig")

    # summary: per source-target-axis
    agg = df.groupby(["source","target","axis"], as_index=False).agg(
        prob=("prob","sum"), min_p=("p_perm","min"), n_pairs=("prob","size"))
    agg.sort_values("prob", ascending=False).to_csv(os.path.join(STEP,"step1_lr_summary_validated.csv"),
                                                   index=False, encoding="utf-8-sig")

    # ---- chord-like plot of the strongest signalling pairs ----
    top = df[(df.p_perm < 0.05) & (df.prob > 0)].head(12)
    if top.empty:
        top = df.head(12)
    nodes = SOURCES + TARGETS
    colors = {"TS":"#4C78A8","MH2":"#F58518","MyoEpi":"#54A24B","M":"#B279A2",
              "PH":"#E45756","GH":"#72B7B2"}
    ang = {n: 2*np.pi*i/len(nodes) for i, n in enumerate(nodes)}
    fig, ax = plt.subplots(figsize=(6.2,6.2))
    for n in nodes:
        th = ang[n]
        ax.add_patch(Wedge((0,0),1.05,np.degrees(th)-8,np.degrees(th)+8,width=0.12,
                           facecolor=colors[n],edgecolor="white",lw=0.5))
        ax.text(1.22*np.cos(th),1.22*np.sin(th),n,ha="center",va="center",fontsize=8)
    mx = top.prob.max()
    for r in top.itertuples():
        a0, a1 = ang[r.source], ang[r.target]
        # quadratic bezier ribbon
        p0 = np.array([np.cos(a0),np.sin(a0)]); p1 = np.array([np.cos(a1),np.sin(a1)])
        ctrl = np.array([0.15*np.cos((a0+a1)/2),0.15*np.sin((a0+a1)/2)])
        path = Path([p0,ctrl,p1], [Path.MOVETO, Path.CURVE3, Path.CURVE3])
        ax.add_patch(PathPatch(path, facecolor="none", edgecolor=colors[r.source],
                               lw=0.6+4.0*r.prob/mx, alpha=0.55))
    ax.set_xlim(-1.35,1.35); ax.set_ylim(-1.35,1.35); ax.axis("off")
    ax.set_title("Niche -> hemocyte LR signalling (top permutation-supported pairs)", fontsize=9)
    fig.savefig(os.path.join(STEP,"step1_chord_validated.png"),dpi=300,bbox_inches="tight")
    fig.savefig(os.path.join(STEP,"step1_chord_validated.pdf"),bbox_inches="tight")
    plt.close(fig)

    # ---- heatmap: source x target, summed probability ----
    piv = df.groupby(["source","target"])["prob"].sum().unstack(fill_value=0).reindex(index=SOURCES,columns=TARGETS,fill_value=0)
    fig, ax = plt.subplots(figsize=(3.6,3.4))
    im = ax.imshow(piv.values, cmap="Reds", aspect="auto")
    ax.set_xticks(range(len(TARGETS))); ax.set_xticklabels(TARGETS)
    ax.set_yticks(range(len(SOURCES))); ax.set_yticklabels(SOURCES)
    for i in range(len(SOURCES)):
        for j in range(len(TARGETS)):
            ax.text(j,i,f"{piv.values[i,j]:.2f}",ha="center",va="center",fontsize=8,
                    color="white" if piv.values[i,j] > piv.values.max()*0.6 else "#333333")
    ax.set_title("Summed LR probability", fontsize=9)
    plt.colorbar(im, ax=ax, fraction=0.045)
    fig.tight_layout()
    fig.savefig(os.path.join(STEP,"step1_lr_heatmap_validated.png"),dpi=300,bbox_inches="tight")
    fig.savefig(os.path.join(STEP,"step1_lr_heatmap_validated.pdf"),bbox_inches="tight")
    plt.close(fig)

    # ---- dot plot of ligand in niche vs receptor in hemocytes ----
    sub = df[df.source.isin(SOURCES) & df.target.isin(TARGETS)].copy()
    sub["pair_label"] = sub["pair"]
    fig, ax = plt.subplots(figsize=(8.0, max(3.5, 0.32*sub.pair_label.nunique())))
    cats = list(dict.fromkeys(sub.pair_label))
    for k, pl in enumerate(cats):
        d = sub[sub.pair_label == pl]
        for _, r in d.iterrows():
            ax.scatter([f"{r.source}->{r.target}"], [k], s=20+300*r.prob/max(df.prob.max(),1e-9),
                       c="#C44E52" if r.p_perm < 0.05 else "#BBBBBB", alpha=0.85,
                       edgecolors="none")
    ax.set_yticks(range(len(cats))); ax.set_yticklabels(cats, fontsize=7)
    ax.tick_params(axis="x", labelrotation=60, labelsize=7)
    ax.set_title("LR communication probability (red = permutation p<0.05)", fontsize=9)
    for s in ax.spines.values(): s.set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(STEP,"step1_lr_dotplot_validated.png"),dpi=300,bbox_inches="tight")
    fig.savefig(os.path.join(STEP,"step1_lr_dotplot_validated.pdf"),bbox_inches="tight")
    plt.close(fig)

    summary = {
        "n_pairs_tested": int(len(df)),
        "n_pairs_p_lt_0.05": int((df.p_perm < 0.05).sum()),
        "top_pairs": df.head(8)[["pair","source","target","prob","p_perm"]].to_dict("records"),
        "note": "CellChat R package unavailable; CellChat-style probability model with label permutation.",
    }
    with open(os.path.join(STEP,"step1_lr_summary_validated.json"),"w",encoding="utf-8") as f:
        json.dump(summary,f,ensure_ascii=False,indent=2)
    print(json.dumps(summary,ensure_ascii=False,indent=2))


if __name__ == "__main__":
    main()
