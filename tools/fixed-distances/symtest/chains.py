"""detect circular equiangular chains (subsets of an orbit of a rotation) inside a configuration, and test the
explicit chain construction S = {0..n-2, n} for n = 6..12."""
import json, os, sys, itertools, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from symmetry import classify, dmat, e_dim
from famcheck import Family
def circ(P, idx):
    A = np.array([[2 * P[i][0], 2 * P[i][1], 1] for i in idx]); b = np.array([P[i][0] ** 2 + P[i][1] ** 2 for i in idx])
    sol, *_ = np.linalg.lstsq(A, b, rcond=None); c = sol[:2]; r = math.sqrt(sol[2] + c @ c)
    res = max(abs(math.dist(P[i], c) - r) for i in idx) / r
    return c, r, res
def chains(P, tol=1e-6):
    """maximal subsets lying on a common circle with angular positions in a common arithmetic progression"""
    n = len(P); found = []
    for m in range(n, 3, -1):
        for idx in itertools.combinations(range(n), m):
            if any(set(idx) <= set(f[0]) for f in found): continue
            c, r, res = circ(P, idx)
            if res > tol: continue
            ang = sorted((math.atan2(P[i][1] - c[1], P[i][0] - c[0]) % (2 * math.pi), i) for i in idx)
            # try every gap as the generator
            gaps = sorted({round(((b[0] - a[0]) % (2 * math.pi)), 9) for a in ang for b in ang if a != b})
            for th in gaps:
                ks = []
                for a in ang:
                    rel = (a[0] - ang[0][0]) % (2 * math.pi)
                    k = next((k for k in range(-40, 41) if abs(((rel - k * th + math.pi) % (2 * math.pi)) - math.pi) < 1e-5), None)
                    if k is None: break
                    ks.append(k)
                if len(ks) == len(ang):
                    S = sorted(ks); S = [k - S[0] for k in S]
                    found.append((idx, S, th)); break
    return found
if __name__ == '__main__':
    D = json.load(open(os.path.join(HERE, '..', 'applet_data.json')))
    for a in ['5,4:473', '5,4:474', '5,4:475', '5,4:477', '5,4:491', '5,4:492', '6,4:23', '6,6:3', '6,6:4', '6,7:0', '6,7:1', '6,7:2', '6,7:3', '6,7:8', '7,7:4']:
        cell, idx = a.split(':'); it = D['cats'][cell]['items'][int(idx)]
        P = [[float(x), float(y)] for x, y in it['coords']]
        print(a, 'chains (points, angular indices S, step/2pi):', [(f[0], f[1], round(f[2] / 2 / math.pi, 5)) for f in chains(P)])
    print()
    print('explicit construction: points at angles k*theta, k in S={0,..,n-2,n}, theta generic')
    for n in range(6, 13):
        S = list(range(n - 1)) + [n]; th = 0.37 + 0.01 * n
        P = [[math.cos(k * th), math.sin(k * th)] for k in S]
        Dm = dmat(P); vals = sorted(Dm[np.triu_indices(n, 1)]); d = 1 + sum(1 for x, y in zip(vals, vals[1:]) if y - x > 1e-9)
        E = list(itertools.combinations(range(n), 2)); uniq = sorted(set(round(v, 9) for v in vals)); cls = [uniq.index(round(Dm[i, j], 9)) for i, j in E]
        F = Family(P, cls); td, sv = F.tangent_dim(F.v)
        print(f'  n={n} d={d} e={e_dim(n, d)} symmetry={classify(P)["type"]} tangent dim={td}')
