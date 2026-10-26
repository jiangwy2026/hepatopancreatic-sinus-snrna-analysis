import io
p = r"D:\桌面\single cell\.venv\Lib\site-packages\scvelo\tools\optimization.py"
src = open(p, encoding="utf8").read()
old = "            gamma[i] = np.linalg.pinv(A.T.dot(A)).dot(A.T.dot(y[:, i]))"
new = "            _g = np.linalg.pinv(A.T.dot(A)).dot(A.T.dot(y[:, i]))\n            gamma[i] = np.asarray(_g).reshape(-1)[0]   # numpy>=2: 1-element array cannot be assigned to a scalar slot"
if "_g = np.linalg.pinv" in src:
    print("already patched")
elif old in src:
    src = src.replace(old, new)
    open(p, "w", encoding="utf8").write(src)
    print("patched", p)
else:
    print("PATTERN NOT FOUND - inspect manually")
