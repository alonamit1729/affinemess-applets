import json, os, sys, itertools, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from symmetry import classify, dmat
from superset import extra_points
from subsets import colouring
D = json.load(open(os.path.join(HERE, '..', 'applet_data.json')))
G = json.load(open(os.path.join(HERE, 'generic_samples.json')))
def nd(Q): return colouring(Q)[1]
for a in sys.argv[1:]:
    cell, idx = a.split(':')
    if a in G: P = G[a][0]['pts']; src = 'generic msolve sample'
    else: it = D['cats'][cell]['items'][int(idx)]; P = [[float(x), float(y)] for x, y in it['coords']]; src = 'catalog sample'
    n = len(P); Dm = dmat(P); vals = sorted(set(round(Dm[i, j], 7) for i, j in itertools.combinations(range(n), 2)))
    m = vals[0]; P = [[x / m, y / m] for x, y in P]; vals = [v / m for v in vals]; d = len(vals)
    ex2 = extra_points(P, vals, new_allowed=2)
    print(f'({cell}) #{idx} [{src}] d={d}: {len(ex2)} candidate extra points with <=2 new distances')
    res = []
    for p in ex2:
        Q = P + [list(p)]; c = classify(Q)['type']
        if c != 'trivial': res.append((1, nd(Q), c, [np.round(p, 5).tolist()]))
    for p, q in itertools.combinations(ex2, 2):
        Q = P + [list(p), list(q)]; k = nd(Q)
        if k <= d + 2:
            c = classify(Q)['type']
            if c != 'trivial': res.append((2, k, c, [np.round(p, 5).tolist(), np.round(q, 5).tolist()]))
    for p, q, r in itertools.combinations(ex2, 3):
        Q = P + [list(p), list(q), list(r)]; k = nd(Q)
        if k <= d + 2:
            c = classify(Q)['type']
            if c != 'trivial': res.append((3, k, c, [np.round(p, 5).tolist(), np.round(q, 5).tolist(), np.round(r, 5).tolist()]))
    res.sort(key=lambda r: (r[1], r[0]))
    for r in res[:8]: print(f'    +{r[0]} point(s) -> {n + r[0]} points, {r[1]} distances, symmetry {r[2]}: {r[3]}')
    if not res: print('    no symmetric superset with <=3 extra points and <=2 new distances')
