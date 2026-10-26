# -*- coding: utf-8 -*-
"""247_step1_communication_sonicparanoid.py

Final Step-1 communication analysis using ONLY orthogroups confirmed by
SonicParanoid 1.3.8 (4 species, sensitive, mmseqs 6.0, vtml80, bitscore 40,
MCL 1.5, in-paralog length diff 0.75).

Orthogroup -> shrimp GeneID mapping comes from step1/new_xp_to_geneid.csv
(new pvannamei proteome XP_ mapped to GCF_003789085.1 GeneIDs).
"""
import os, json
import numpy as np, pandas as pd, anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, PathPatch
from matplotlib.path import Path

OUT = r"D:\桌面\single cell\analysis_output"
STEP = os.path.join(OUT, "step1")
os.environ.setdefault("MPLCONFIGDIR", os.path.join(OUT, ".mplconfig"))
RNG = np.random.default_rng(42)

SOURCES = ["TS", "MH2", "MyoEpi", "M"]
TARGETS = ["PH", "GH"]
ALL = SOURCES + TARGETS

# (ligand name, ligand GeneIDs, receptor name, receptor GeneIDs, axis)
LR = [
 ("JAG1", ["113821858"], "NOTCH1", ["113802681"], "Notch"),
 ("DLL1", ["113815864"], "NOTCH1", ["113802681"], "Notch"),
 ("WNT5A", ["113805832","113811410","113806005","113827257"], "FZD4", ["113819177","113826275"], "Wnt"),
 ("WNT5A", ["113805832","113811410","113806005","113827257"], "FZD5", ["113826277"], "Wnt"),
 ("THBS4", ["113823118"], "ITGAV", ["113800203","113821319"], "ECM-integrin"),
 ("THBS4", ["113823118"], "ITGB1", ["113830322","113805715"], "ECM-integrin"),
 ("LAMA1", ["113823054","113828881"], "ITGAV", ["113800203","113821319"], "ECM-integrin"),
 ("LAMA1", ["113823054","113828881"], "ITGB1", ["113830322","113805715"], "ECM-integrin"),
 ("LAMA1", ["113823054","113828881"], "ITGA4", ["113827435"], "ECM-integrin"),
 ("POSTN", ["113817567","113805211"], "ITGAV", ["113800203","113821319"], "ECM-integrin"),
 ("POSTN", ["113817567","113805211"], "ITGB1", ["113830322","113805715"], "ECM-integrin"),
 ("VWF", ["113824171"], "ITGAV", ["113800203","113821319"], "ECM-integrin"),
 ("VWF", ["113824171"], "ITGB1", ["113830322","113805715"], "ECM-integrin"),
]


def main():
    a = ad.read_h5ad(os.path.join(STEP, "step1_niche.h5ad"))
    a.X = a.layers["lognorm"].copy()
    g2v = {str(g).replace("ncbi_", ""): v for g, v in zip(a.var["gid"], a.var_names)}
    labels = a.obs["cell_type"].astype(str).values

    def group_expr(gids):
        idx = [g2v[g] for g in gids if g in g2v]
        if not idx:
            return np.zeros(a.n_obs)
        X = np.asarray(a[:, idx].X.todense(), dtype=np.float64).mean(axis=1)
        return X

    rows = []
    n_perm = 1000
    for lig, lgids, rec, rgids, axis in LR:
        xl, xr = group_expr(lgids), group_expr(rgids)
        ml = {g: float(xl[labels == g].mean()) for g in ALL}
        mr = {g: float(xr[labels == g].mean()) for g in ALL}
        for s in SOURCES:
            for t in TARGETS:
                prob = ml[s] * mr[t]
                null = np.empty(n_perm)
                for k in range(n_perm):
                    perm = RNG.permutation(labels)
                    null[k] = xl[perm == s].mean() * xr[perm == t].mean()
                p = (np.sum(null >= prob) + 1) / (n_perm + 1)
                rows.append({"ligand": lig, "receptor": rec, "axis": axis, "pair": f"{lig} -> {rec}",
                             "source": s, "target": t, "ligand_mean": ml[s], "receptor_mean": mr[t],
                             "prob": prob, "p_perm": p, "ligand_genes": ";".join(lgids),
                             "receptor_genes": ";".join(rgids)})
    df = pd.DataFrame(rows).sort_values("prob", ascending=False)
    df.to_csv(os.path.join(STEP, "step1_lr_communication_sonicparanoid.csv"), index=False, encoding="utf-8-sig")

    agg = df.groupby(["source","target","axis"], as_index=False).agg(prob=("prob","sum"), min_p=("p_perm","min"))
    agg.sort_values("prob", ascending=False).to_csv(os.path.join(STEP,"step1_lr_summary_sonicparanoid.csv"),
                                                   index=False, encoding="utf-8-sig")

    top = df[(df.p_perm < 0.05) & (df.prob > 0)].head(12)
    if top.empty: top = df.head(12)
    nodes = SOURCES + TARGETS
    colors = {"TS":"#4C78A8","MH2":"#F58518","MyoEpi":"#54A24B","M":"#B279A2","PH":"#E45756","GH":"#72B7B2"}
    ang = {n: 2*np.pi*i/len(nodes) for i, n in enumerate(nodes)}
    fig, ax = plt.subplots(figsize=(6.2,6.2))
    for n in nodes:
        th = ang[n]
        ax.add_patch(Wedge((0,0),1.05,np.degrees(th)-8,np.degrees(th)+8,width=0.12,
                           facecolor=colors[n], edgecolor="white", lw=0.5))
        ax.text(1.22*np.cos(th),1.22*np.sin(th),n,ha="center",va="center",fontsize=8)
    mx = max(top.prob.max(), 1e-9)
    for r in top.itertuples():
        a0, a1 = ang[r.source], ang[r.target]
        p0 = np.array([np.cos(a0),np.sin(a0)]); p1 = np.array([np.cos(a1),np.sin(a1)])
        ctrl = np.array([0.15*np.cos((a0+a1)/2),0.15*np.sin((a0+a1)/2)])
        path = Path([p0,ctrl,p1],[Path.MOVETO,Path.CURVE3,Path.CURVE3])
        ax.add_patch(PathPatch(path, facecolor="none", edgecolor=colors[r.source],
                               lw=0.6+4.0*r.prob/mx, alpha=0.55))
    ax.set_xlim(-1.35,1.35); ax.set_ylim(-1.35,1.35); ax.axis("off")
    ax.set_title("Niche -> hemocyte signalling (SonicParanoid-confirmed orthogroups)", fontsize=9)
    fig.savefig(os.path.join(STEP,"step1_chord_sonicparanoid.png"),dpi=300,bbox_inches="tight")
    fig.savefig(os.path.join(STEP,"step1_chord_sonicparanoid.pdf"),bbox_inches="tight")
    plt.close(fig)

    piv = df.groupby(["source","target"])["prob"].sum().unstack(fill_value=0).reindex(index=SOURCES,columns=TARGETS,fill_value=0)
    fig, ax = plt.subplots(figsize=(3.6,3.4))
    im = ax.imshow(piv.values, cmap="Reds", aspect="auto")
    ax.set_xticks(range(len(TARGETS))); ax.set_xticklabels(TARGETS)
    ax.set_yticks(range(len(SOURCES))); ax.set_yticklabels(SOURCES)
    for i in range(len(SOURCES)):
        for j in range(len(TARGETS)):
            ax.text(j,i,f"{piv.values[i,j]:.3f}",ha="center",va="center",fontsize=8,
                    color="white" if piv.values[i,j] > piv.values.max()*0.6 else "#333333")
    ax.set_title("Summed LR probability (SonicParanoid axes)", fontsize=9)
    plt.colorbar(im, ax=ax, fraction=0.045)
    fig.tight_layout(); fig.savefig(os.path.join(STEP,"step1_lr_heatmap_sonicparanoid.png"),dpi=300,bbox_inches="tight")
    fig.savefig(os.path.join(STEP,"step1_lr_heatmap_sonicparanoid.pdf"),bbox_inches="tight"); plt.close(fig)

    summary = {"n_tests": int(len(df)), "n_p_lt_0.05": int((df.p_perm < 0.05).sum()),
               "axes": df.groupby("axis").prob.sum().sort_values(ascending=False).to_dict(),
               "top": df.head(10)[["pair","source","target","prob","p_perm"]].to_dict("records")}
    with open(os.path.join(STEP,"step1_lr_summary_sonicparanoid.json"),"w",encoding="utf-8") as f:
        json.dump(summary,f,ensure_ascii=False,indent=2)
    print(json.dumps(summary,ensure_ascii=False,indent=2))


if __name__ == "__main__":
    main()
