"""Collect a finished run of run_level.py.

    python3 finish_level.py N OUTDIR CANDIDATES.json TARGET_DIR

Writes TARGET_DIR/levelN.jsonl.gz (all certificates, in candidate order; the chunk files are gzip
members, so the file is their concatenation) and TARGET_DIR/LN.json (the candidates that are kept),
and prints the counts.
"""
import glob, hashlib, json, os, sys
from collections import Counter

n, outdir, source, target = int(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4]
total = len(json.load(open(source))["patterns"])
names = sorted(glob.glob(os.path.join(outdir, "c[0-9][0-9][0-9][0-9].json")))
summaries = [json.load(open(name)) for name in names]
assert [s["chunk"] for s in summaries] == list(range(len(summaries))), "chunks missing"
assert sum(s["candidates"] for s in summaries) == total, "run not finished"
kept = [p for s in summaries for p in s["kept"]]
degrees = Counter()
for s in summaries:
    degrees.update({int(d): c for d, c in s["degrees"].items()})
certificates = sum(s["certificates"] for s in summaries)
os.makedirs(target, exist_ok=True)
sha = hashlib.sha256()
with open(os.path.join(target, f"level{n}.jsonl.gz"), "wb") as out:
    for name in names:
        with open(name[:-5] + ".jsonl.gz", "rb") as fh:
            data = fh.read()
        out.write(data)
        sha.update(data)
with open(os.path.join(target, f"L{n}.json"), "w") as fh:
    json.dump({"n": n, "count": len(kept), "patterns": kept}, fh, separators=(",", ":"))
tally = lambda ps: dict(sorted(Counter(max(p) + 1 for p in ps).items()))
print(f"{total} candidates on {n} points: {certificates} certificates, {len(kept)} kept")
print("kept, by number of colours:", tally(kept))
print("degree of the product -> certificates:", dict(sorted(degrees.items())))
print(f"certificates: {sum(s['bytes'] for s in summaries):,} bytes as JSON lines,",
      f"{os.path.getsize(os.path.join(target, f'level{n}.jsonl.gz')):,} bytes compressed")
print(f"search and check: {sum(s['seconds'] for s in summaries) / 3600:.1f} core-hours")
print(f"sha256 of level{n}.jsonl.gz: {sha.hexdigest()}")
