"""Certificates for a whole level, in chunks, resumable.

    python3 run_level.py N CANDIDATES.json OUTDIR [WORKERS]

CANDIDATES.json holds {"patterns": [...]} (all candidates on N points, in a fixed order).  For chunk i
of 4,000 candidates the run writes OUTDIR/cNNNN.jsonl.gz (the certificates, in candidate order) and
OUTDIR/cNNNN.json (the candidates without a certificate, and counts).  Chunks that exist are skipped.
Every certificate is checked with check.py before it is written.
"""
import gzip, json, os, sys, time
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import check
import fast

CHUNK = 4000


def work(job):
    n, index, patterns, outdir = job
    t0 = time.time()
    lines, kept, degrees = [], [], {}
    for p in patterns:
        cert, d = fast.search(n, p)
        if cert is None:
            kept.append(p)
        else:
            check.check(cert)
            lines.append(json.dumps(cert, separators=(",", ":")))
            degrees[d] = degrees.get(d, 0) + 1
    data = ("\n".join(lines) + "\n").encode() if lines else b""
    name = os.path.join(outdir, f"c{index:04d}")
    with open(name + ".tmp", "wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, filename="") as fh:
            fh.write(data)
    os.replace(name + ".tmp", name + ".jsonl.gz")
    summary = {"chunk": index, "candidates": len(patterns), "certificates": len(lines), "kept": kept,
               "degrees": degrees, "bytes": len(data), "seconds": round(time.time() - t0, 1)}
    with open(name + ".json.tmp", "w") as fh:
        json.dump(summary, fh, separators=(",", ":"))
    os.replace(name + ".json.tmp", name + ".json")
    return index, len(patterns), len(lines), len(kept), summary["seconds"]


if __name__ == "__main__":
    n, source, outdir = int(sys.argv[1]), sys.argv[2], sys.argv[3]
    workers = int(sys.argv[4]) if len(sys.argv) > 4 else 2
    os.makedirs(outdir, exist_ok=True)
    patterns = json.load(open(source))["patterns"]
    jobs = [(n, i // CHUNK, patterns[i:i + CHUNK], outdir) for i in range(0, len(patterns), CHUNK)
            if not os.path.exists(os.path.join(outdir, f"c{i // CHUNK:04d}.json"))]
    print(f"{len(patterns)} candidates, {len(jobs)} chunks to do", flush=True)
    t0, done = time.time(), 0
    with Pool(workers) as pool:
        for index, m, c, k, s in pool.imap_unordered(work, jobs):
            done += 1
            print(f"chunk {index}: {c} certificates, {k} kept, {s:.0f} s   [{done}/{len(jobs)} chunks, {time.time() - t0:.0f} s]", flush=True)
    print("finished", flush=True)
