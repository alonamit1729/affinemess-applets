"""For a counterexample sample: find extra points p whose distances to all n points lie in the configuration's
distance set (possibly also allowing one new value), and the symmetry of the enlarged configurations."""
import json, os, sys, itertools, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from symmetry import classify, dmat
def circle_pts(a, ra, b, rb):
    d = math.dist(a, b)
    if d < 1e-12 or d > ra + rb + 1e-9 or d < abs(ra - rb) - 1e-9: return []
    x = (ra * ra - rb * rb + d * d) / (2 * d); h2 = ra * ra - x * x
    h = math.sqrt(max(h2, 0)); ux, uy = (b[0] - a[0]) / d, (b[1] - a[1]) / d
    px, py = a[0] + x * ux, a[1] + x * uy
    return [(px - h * uy, py + h * ux), (px + h * uy, py - h * ux)]
def extra_points(P, vals, tol=1e-6, new_allowed=0):
    n = len(P); sc = max(vals); out = []
    for i, j in itertools.combinations(range(n), 2):
        for ra in vals:
            for rb in vals:
                for p in circle_pts(P[i], ra, P[j], rb):
                    if any(math.dist(p, q) < tol * sc for q in P): continue
                    ds = [math.dist(p, q) for q in P]
                    bad = [x for x in ds if min(abs(x - v) for v in vals) > tol * sc]
                    if len(bad) <= new_allowed and not any(math.dist(p, q) < tol * sc for q in out):
                        out.append(p)
    return out
if __name__ == '__main__':
    D = json.load(open(os.path.join(HERE, '..', 'applet_data.json')))
    for a in sys.argv[1:]:
        cell, idx = a.split(':'); it = D['cats'][cell]['items'][int(idx)]
        P = [[float(x), float(y)] for x, y in it['coords']]; n = len(P); Dm = dmat(P)
        vals = sorted(set(round(Dm[i, j], 7) for i, j in itertools.combinations(range(n), 2)))
        m = vals[0]; P = [[x / m, y / m] for x, y in P]; vals = [v / m for v in vals]
        ex = extra_points(P, vals)
        print(f'({cell}) #{idx}: d={len(vals)}, extra points with all distances in the distance set: {len(ex)}')
        for p in ex:
            Q = P + [list(p)]; c = classify(Q)
            print(f'    p={np.round(p, 5).tolist()}  symmetry of the {n+1}-point set: {c["type"]}')
        if len(ex) >= 2:
            Q = P + [list(p) for p in ex]; c = classify(Q); Dq = dmat(Q)
            dq = len(set(round(Dq[i, j], 6) for i, j in itertools.combinations(range(len(Q)), 2)))
            print(f'    all extra points together: {len(Q)} points, {dq} distances, symmetry {c["type"]}')
