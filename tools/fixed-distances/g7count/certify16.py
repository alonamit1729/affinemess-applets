"""Certification of each 16-point set found: msolve real root isolation on its 16-point colour pattern (coordinate
system of realize.py: p0=(0,0), p1=(x1,0), colour of edge (0,1) has value 1, free values s1..s6, Rabinowitsch t),
precision 120.  Output: certify16.json with the number of real solutions and the similarity classes among them."""
import json, sys, time, os
from msolve_util import *
from extend import *
R = json.load(open(sys.argv[1]))['sets']
out = []
for i, r in enumerate(R):
    col = tuple(r['pattern'])
    t0 = time.time()
    V, sols, dt = realize_pattern(16, col, precision=120, timeout=20000)
    rec = {'set': i, 'k': max(col) + 1, 'time': round(dt, 1)}
    if isinstance(sols, str):
        rec['status'] = sols
    else:
        rec['nsol'] = len(sols)
        C = Classes()
        from fractions import Fraction
        for sol in sols:
            vv = {}
            for name, (a, b) in zip(V, sol):
                m = (a + b) / 2
                vv[name] = mp.mpf(m.numerator) / mp.mpf(m.denominator)
            P = [(mp.mpf(0), mp.mpf(0)), (vv['x1'], mp.mpf(0))] + [(vv[f'x{j}'], vv[f'y{j}']) for j in range(2, 16)]
            assert pattern(P) == col
            C.add(P)
        rec['nclasses'] = len(C.reps)
        orig = [(mp.mpf(x), mp.mpf(y)) for x, y in r['pts']]
        rec['matches_found_set'] = [congruent(c[0], orig) for c in C.reps]
        rec['max_width'] = float(max(b - a for sol in sols for a, b in sol))
    print(rec, flush=True); out.append(rec)
    json.dump(out, open('certify16.json', 'w'), indent=1)
