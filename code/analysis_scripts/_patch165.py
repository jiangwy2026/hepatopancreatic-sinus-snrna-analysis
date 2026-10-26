p = r"D:\桌面\single cell\analysis_output\scripts\165_final_annotation_figs.py"
src = open(p, encoding="utf8").read().replace("\r\n", "\n")
old = '''    P = pd.read_csv(os.path.join(OUT, "C0_metabolic_programs.csv"), index_col=0)
    P = P.rename(columns=NEW)
    P = P[[c for c in ORDER if c in P.columns]]
    P.index.name = "program"
    P.to_csv(os.path.join(OUT, "metabolic_programs_by_cluster.csv"))'''
new = '''    P = pd.read_csv(os.path.join(OUT, "C0_metabolic_programs.csv"), index_col=0)
    PMAP = {"MH1(C0)": "TS", "GH": "GH", "PH": "PH", "MH2": "MH2", "Epi": "PM-like",
            "F": "F", "R": "R", "B": "B", "Neu8": "MyoEpi", "Neu9": "M"}
    P = P.rename(columns=PMAP)
    P = P[["PH", "GH", "TS", "MH2", "F", "R", "B", "M", "MyoEpi", "PM-like"]]
    P.index.name = "program"
    P.to_csv(os.path.join(OUT, "metabolic_programs_by_cluster.csv"))
    # 热图（行内 z-score）
    Z = P.sub(P.mean(1), axis=0).div(P.std(1).replace(0, np.nan), axis=0)
    fig, ax = plt.subplots(figsize=(5.6, 2.2))
    im = ax.imshow(Z.values, cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
    ax.set_xticks(range(P.shape[1])); ax.set_xticklabels(P.columns, rotation=45, ha="right", fontsize=6.5)
    ax.set_yticks(range(P.shape[0])); ax.set_yticklabels(P.index, fontsize=6.2)
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            if abs(Z.values[i, j]) > 1.0:
                ax.text(j, i, "%.1f" % Z.values[i, j], fontsize=5, ha="center", va="center",
                        color="white" if abs(Z.values[i, j]) > 1.6 else "#333")
    ax.set_title("Metabolic / biosynthetic programs by cluster (row z-score)", fontsize=7.5)
    cb = fig.colorbar(im, ax=ax, fraction=0.03, shrink=0.7, pad=0.02)
    cb.ax.tick_params(labelsize=4.8); cb.outline.set_visible(False)
    for s in ax.spines.values(): s.set_visible(False)
    fig.savefig(os.path.join(FIG, "Fig_metabolic_programs_preview.png"), dpi=400, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)'''
assert old in src; src = src.replace(old, new, 1)
open(p, "w", encoding="utf8", newline="\n").write(src); print("patched")
