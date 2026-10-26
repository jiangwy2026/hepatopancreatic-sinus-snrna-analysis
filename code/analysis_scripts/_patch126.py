p = r"D:\桌面\single cell\analysis_output\scripts\126_dropped_clusters.py"
src = open(p, encoding="utf8").read().replace("\r\n", "\n")
src = src.replace("sigK[g] = np.asarray(K[:, gi].mean(1)).ravel()", "sigK[g] = np.asarray(K[:, gi].X.mean(1)).ravel()")
src = src.replace('S[g] = np.asarray(D2[:, gi].mean(1)).ravel()', 'S[g] = np.asarray(D2[:, gi].X.mean(1)).ravel()')
src = src.replace('S[name] = np.asarray(D2[:, gi].mean(1)).ravel()', 'S[name] = np.asarray(D2[:, gi].X.mean(1)).ravel()')
open(p, "w", encoding="utf8", newline="\n").write(src); print("patched")
