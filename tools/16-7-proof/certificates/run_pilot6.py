"""Six points: the pilot.  A random sample of 1,060 of the 537,230 six-point candidates
(pilot6_sample.json: index in the sorted list of candidates, and the pattern).

    python3 run_pilot6.py          (about one minute on two cores)

Writes pilot6.jsonl.gz (the certificates found) and pilot6_kept.json (the candidates without one).
"""
import gzip, json, os, sys, time
from collections import Counter
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import check
from find import search


def work(item):
    index, pattern = item
    t0 = time.time()
    cert, hf, degree = search(6, tuple(pattern), dmax=6)
    t1 = time.time()
    if cert is not None:
        check.check(cert)
    return index, pattern, cert, degree, t1 - t0, time.time() - t1


if __name__ == "__main__":
    sample = json.load(open(os.path.join(HERE, "pilot6_sample.json")))["candidates"]
    t0 = time.time()
    with Pool(2) as pool:
        results = sorted(pool.imap_unordered(work, sample, chunksize=4))
    certs = [r for r in results if r[2] is not None]
    kept = [r for r in results if r[2] is None]
    sizes = sorted(len(json.dumps(r[2], separators=(",", ":"))) for r in certs)
    print(f"{len(results)} candidates: {len(certs)} with a certificate, {len(kept)} kept   ({time.time() - t0:.0f} s)")
    print("degree of the product -> certificates:", dict(sorted(Counter(r[3] for r in certs).items())))
    print(f"bytes per certificate: median {sizes[len(sizes) // 2]}, mean {sum(sizes) // len(sizes)}, largest {sizes[-1]}")
    print(f"time per candidate: search {1000 * sum(r[4] for r in results) / len(results):.0f} ms,",
          f"check {1000 * sum(r[5] for r in certs) / len(certs):.1f} ms per certificate")
    data = "".join(json.dumps(r[2], separators=(",", ":")) + "\n" for r in certs).encode()
    with open(os.path.join(HERE, "pilot6.jsonl.gz"), "wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, filename="") as fh:
            fh.write(data)
    with open(os.path.join(HERE, "pilot6_kept.json"), "w") as fh:
        json.dump({"n": 6, "count": len(kept), "patterns": [r[1] for r in kept]}, fh, separators=(",", ":"))
