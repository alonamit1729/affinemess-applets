"""Five points: certificates for all discarded candidates.  Writes level5.jsonl and L5.json.

    python3 run_level5.py          (about 8 minutes; almost all of it is spent on the patterns that are kept)
"""
import json, os, sys, time
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))          # fewdist.py is one folder up
sys.path.insert(0, HERE)
import fewdist as fd
import check
from find import search

four = fd.all_patterns(4)
L4 = [p for p in four if p != (0,) * 6]            # the one-coloured K4 is not realisable (certificate below)
raw5 = sorted(fd.add_a_point(L4, allowed=L4))
t0 = time.time()
certs, kept = [], []
cert4, _, _ = search(4, (0,) * 6)
check.check(cert4)
for p in raw5:
    cert, _, _ = search(5, p, dmax=6)
    if cert is None:
        kept.append(p)
    else:
        check.check(cert)
        certs.append(cert)
tally = lambda ps: dict(sorted(Counter(max(p) + 1 for p in ps).items()))
print(f"{len(raw5)} candidates; {len(certs)} with a certificate, by colours {tally([c['pattern'] for c in certs])};")
print(f"{len(kept)} kept, by colours {tally(kept)}   ({time.time() - t0:.0f} s)")
with open(os.path.join(HERE, "level5.jsonl"), "w") as fh:
    for c in [cert4] + certs:
        fh.write(json.dumps(c, separators=(",", ":")) + "\n")
with open(os.path.join(HERE, "L5.json"), "w") as fh:
    json.dump({"n": 5, "count": len(kept), "patterns": [list(p) for p in kept]}, fh, separators=(",", ":"))
