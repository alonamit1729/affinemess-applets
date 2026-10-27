"""verify that the symmetric supersets found by superset2 are genuine families: tangent dimension of the enlarged
colouring system and a numerical trace keeping the symmetry; then msolve confirmation of the enlarged system."""
import json, os, sys, itertools, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from symmetry import classify, dmat, e_dim
from famcheck import Family
from subsets import colouring
from superset import extra_points
D = json.load(open(os.path.join(HERE, '..', 'applet_data.json')))
G = json.load(open(os.path.join(HERE, 'generic_samples.json')))
cases = {'5,4:474': [[2.3615, -2.18766]], '5,4:475': [[1.4191, -0.53987]], '5,4:477': [[-2.34733, -3.0657]], '6,4:23': [[1.60332, 0.95507]],
         '6,6:4': [[-0.17531, 1.65841], [-1.04134, 1.15841]], '6,7:1': [[0.80714, -1.64714]], '6,7:2': [[0.79414, 0.80417]],
         '6,7:3': [[-0.0542, -0.99853], [1.75147, 1.59916]], '6,7:8': [[0.23827, 0.9712], [1.44874, 1.0226], [2.82044, -0.69195]]}
out = {}
for a, extra in cases.items():
    cell, idx = a.split(':'); it = D['cats'][cell]['items'][int(idx)]
    P = G[a][0]['pts'] if a in G else [[float(x), float(y)] for x, y in it['coords']]; n = len(P); Dm = dmat(P)
    m = min(Dm[i, j] for i, j in itertools.combinations(range(n), 2)); P = [[x / m, y / m] for x, y in P]
    vals = sorted(set(round(Dm[i, j] / m, 7) for i, j in itertools.combinations(range(n), 2)))
    cand = extra_points(P, vals, new_allowed=2)
    extra = [min(cand, key=lambda q: math.dist(q, e)) for e in extra]
    Q = np.array(P + [list(e) for e in extra]); cls, dq = colouring(Q)
    F = Family(Q, cls); v, ok = F.project(F.v); td, sv = F.tangent_dim(v)
    T = F.tangents(v); syms = set(); ends = []
    for u in T[:1]:
        for sgn in (1, -1):
            tr, vend = F.trace(v, sgn * u, steps=200)
            syms |= {t[1] for t in tr if t[1] not in ('boundary', 'lost', 'singular')}; ends.append(classify(F.points(vend))['type'])
    nq = len(Q)
    print(f'{a} + {len(extra)} point(s): {nq} points, {dq} distances, e={e_dim(nq, dq)}, residual {abs(F.res(v)).max():.1e}, tangent dim {td}, symmetry along trace {syms}, ends {ends}')
    out[a] = {'points': F.points(v).tolist(), 'cls': cls, 'n': nq, 'd': dq, 'e': e_dim(nq, dq), 'tangent_dim': td, 'syms': sorted(syms)}
json.dump(out, open(os.path.join(HERE, 'superset_families.json'), 'w'), indent=1)
