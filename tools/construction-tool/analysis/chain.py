"""Practical constructibility: can a player build the configuration from a base
(two points / a lattice / a regular polygon with its centre) by repeatedly adding
a configuration point that lies on two distinct objects already available:
  - a circle about an available point with radius equal to a distance already
    present among available points (compass),
  - a line through two available points,
  - the perpendicular bisector of two available points?
Usage: python3 chain.py rigid.json out.json
"""
import json, sys, math, cmath, itertools
from multiprocessing import Pool

TOL = 1e-8


FLAGS = dict(bis=True, lines=True, kmax=None, centre_vertex=False)

def close(S, cfg, use_bisector=None):
    use_bisector = FLAGS['bis'] if use_bisector is None else use_bisector
    """S: list of complex (available points, may include auxiliary). cfg: list of complex targets.
    Returns (set of reached cfg indices, number of rounds)."""
    S = list(S)
    reached = set(i for i, p in enumerate(cfg) if any(abs(p - q) < TOL for q in S))
    rounds = 0
    while True:
        D = [abs(a - b) for a, b in itertools.combinations(S, 2)]
        Dset = D
        added = []
        for i, p in enumerate(cfg):
            if i in reached:
                continue
            objs = []
            def addobj(o):
                for q in objs:
                    if q[0] != o[0]:
                        continue
                    if o[0] == 'c':
                        if abs(q[1] - o[1]) < 1e-9 and abs(q[2] - o[2]) < 1e-6:
                            return
                    else:
                        # lines: unit direction u (mod sign) and offset
                        d = abs((q[1].conjugate() * o[1]).imag)
                        if d < 1e-6:
                            # same direction; same line if o's point lies on q
                            if abs((q[1].conjugate() * (o[2] - q[2])).imag) < 1e-6:
                                return
                objs.append(o)
            for k, a in enumerate(S):
                r = abs(p - a)
                if r < TOL:
                    continue
                if any(abs(r - d) < TOL for d in Dset):
                    addobj(('c', a, r))
            for (k, a), (l, b) in itertools.combinations(enumerate(S), 2):
                ab = b - a
                if abs(ab) < TOL:
                    continue
                u = ab / abs(ab)
                if FLAGS['lines'] and abs((u.conjugate() * (p - a)).imag) < TOL:
                    addobj(('l', u, a))
                if use_bisector and abs(abs(p - a) - abs(p - b)) < TOL:
                    addobj(('l', u * 1j, (a + b) / 2))
            if len(objs) >= 2:
                added.append(i)
        if not added:
            return reached, rounds
        rounds += 1
        for i in added:
            reached.add(i)
            S.append(cfg[i])


def polygon(a, b, m, k, sign):
    """regular m-gon in which a, b are vertices k steps apart (orientation sign). returns centre, vertices"""
    th = 2 * math.pi * k / m * sign
    # centre c with (b - c) = e^{i th} (a - c)
    w = cmath.exp(1j * th)
    c = (b - w * a) / (1 - w)
    v = [c + (a - c) * cmath.exp(1j * 2 * math.pi * j / m) for j in range(m)]
    return c, v


def lattice_hits(a, b, cfg, kind):
    u = b - a
    v = u * (cmath.exp(1j * math.pi / 3) if kind == 'tri' else 1j)
    hits = []
    for p in cfg:
        z = p - a
        # solve z = s u + t v
        det = (u.conjugate() * v).imag
        s = (z.conjugate() * v).imag / det * -1
        t = (u.conjugate() * z).imag / det
        s = ((z * v.conjugate()).imag) / ((u * v.conjugate()).imag)
        t = ((z * u.conjugate()).imag) / ((v * u.conjugate()).imag)
        if abs(s - round(s)) < 1e-6 and abs(t - round(t)) < 1e-6:
            hits.append(p)
    return hits


def analyse(o):
    cfg = [complex(float(x), float(y)) for x, y in o['coords']]
    n = len(cfg)
    pairs = list(itertools.combinations(range(n), 2))
    lat = []
    for kind in ('tri', 'sq'):
        for i, j in itertools.permutations(range(n), 2):
            h = lattice_hits(cfg[i], cfg[j], cfg, kind)
            if len(h) == n:
                return dict(how='lattice-' + kind, rounds=0)
            if len(h) >= 2:
                lat.append((kind, h))
    MS = list(range(3, 19)) + [20, 24, 30]
    for m in MS:
        for i, j in pairs:
            for k in range(1, (FLAGS['kmax'] or m // 2) + 1):
                for sgn in (1, -1):
                    c, V = polygon(cfg[i], cfg[j], m, k, sgn)
                    base = V + [c]
                    if all(any(abs(p - q) < TOL for q in base) for p in cfg):
                        return dict(how='polygon-subset', m=m, rounds=0)
    best = None
    for kind, h in lat:
        r, rd = close(h, cfg)
        if len(r) == n and (best is None or rd < best[1]):
            best = (kind, rd)
    if best:
        return dict(how='lattice-' + best[0] + '+chain', rounds=best[1])
    tbest = None
    for i, j in pairs:
        r, rd = close([cfg[i], cfg[j]], cfg)
        if len(r) == n and (tbest is None or rd < tbest):
            tbest = rd
    if tbest is not None:
        return dict(how='pair+chain', rounds=tbest)
    for m in MS:
        pb = None
        for i, j in pairs:
            for k in range(1, (FLAGS['kmax'] or m // 2) + 1):
                for sgn in (1, -1):
                    c, V = polygon(cfg[i], cfg[j], m, k, sgn)
                    r, rd = close(V + [c], cfg)
                    if len(r) == n and (pb is None or rd < pb):
                        pb = rd
        if pb is not None:
            return dict(how='polygon+chain', m=m, rounds=pb)
    return dict(how=None)


def work(o):
    r = analyse(o)
    r.update(cell=o['cell'], i=o['i'], label=o['label'], frame=o['frame'], miss=o['miss'])
    return r


if __name__ == '__main__':
    R = json.load(open(sys.argv[1]))
    with Pool(2) as P:
        out = []
        for k, r in enumerate(P.imap_unordered(work, R, chunksize=4)):
            out.append(r)
            if k % 100 == 0:
                print(k, len(R), flush=True)
                json.dump(out, open(sys.argv[2], 'w'))
    json.dump(out, open(sys.argv[2], 'w'))
    print('done', len(out))
