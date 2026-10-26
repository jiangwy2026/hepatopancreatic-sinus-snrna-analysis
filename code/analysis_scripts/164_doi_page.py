import io, sys, re, requests
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
h = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}
for url in ("https://doi.org/10.1016/j.fsi.2026.111494", "https://www.sciencedirect.com/science/article/pii/S1050464826001117"):
    try:
        r = requests.get(url, headers=h, timeout=30, allow_redirects=True)
        print(url, "->", r.status_code, "| final:", r.url[:100], "| len:", len(r.text))
        links = set(re.findall(r'(https?://[^"\']*?(?:mmc\d+|supplementary|appendix)[^"\']*)', r.text, re.I))
        print("   补充材料链接:", list(links)[:5] if links else "(无)")
    except Exception as e:
        print(url, "失败:", repr(e)[:100])
