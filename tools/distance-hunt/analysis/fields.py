"""For each rigid configuration: refine to high precision, then find the minimal
polynomial of every squared-distance ratio and probe its Galois group by
factoring modulo many primes (Frobenius cycle types).

Usage: python3 fields.py rigid.json out.json [only_missing]
"""
import json, sys, math
from multiprocessing import Pool
import mpmath as mp
import flint

DPS = 1300
PRIMES = [p for p in range(3, 60000) if flint.fmpz(p).is_prime()][:1500]


def refine(coords, cls, n):
    mp.mp.dps = 30
    z = [mp.mpc(mp.mpf(x), mp.mpf(y)) for x, y in coords]
    z0, z1 = z[0], z[1]
    w = [(q - z0) / (z1 - z0) for q in z]
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    c01 = cls[0]
    classes = sorted(set(cls))
    free = [c for c in classes if c != c01]
    # unknown vector: x2,y2,...,x_{n-1},y_{n-1}, t_c for c in free
    v = []
    for q in w[2:]:
        v += [q.real, q.imag]
    t0 = {}
    for (i, j), c in zip(pairs, cls):
        t0.setdefault(c, abs(w[i] - w[j]) ** 2)
    v += [t0[c] for c in free]
    m = len(v)
    dps = 30
    for it in range(40):
        mp.mp.dps = dps + 20
        v = [mp.mpf(a) for a in v]
        P = [(mp.mpf(0), mp.mpf(0)), (mp.mpf(1), mp.mpf(0))] + [(v[2 * k], v[2 * k + 1]) for k in range(n - 2)]
        T = {c01: mp.mpf(1)}
        for idx, c in enumerate(free):
            T[c] = v[2 * (n - 2) + idx]
        F = []
        J = []
        for (i, j), c in zip(pairs, cls):
            dx = P[i][0] - P[j][0]; dy = P[i][1] - P[j][1]
            F.append(dx * dx + dy * dy - T[c])
            row = [mp.mpf(0)] * m
            if i >= 2:
                row[2 * (i - 2)] += 2 * dx; row[2 * (i - 2) + 1] += 2 * dy
            if j >= 2:
                row[2 * (j - 2)] -= 2 * dx; row[2 * (j - 2) + 1] -= 2 * dy
            if c != c01:
                row[2 * (n - 2) + free.index(c)] -= 1
            J.append(row)
        Jm = mp.matrix(J); Fm = mp.matrix(F)
        A = Jm.T * Jm; b = Jm.T * Fm
        delta = mp.lu_solve(A, b)
        v = [v[k] - delta[k] for k in range(m)]
        res = max(abs(f) for f in F)
        step = max(abs(delta[k]) for k in range(m))
        if step < mp.mpf(10) ** (-(DPS + 10)):
            break
        if step < mp.mpf(10) ** (-dps // 2):
            dps = min(2 * dps, DPS + 40)
    mp.mp.dps = DPS
    T = {c01: mp.mpf(1)}
    for idx, c in enumerate(free):
        T[c] = v[2 * (n - 2) + idx]
    # rank check (rigidity): smallest singular value proxy
    return T, float(res), float(step)


def minpoly(x, dmax):
    mp.mp.dps = DPS
    prec_digits = DPS - 100
    C = mp.mpf(10) ** prec_digits
    rows = []
    pw = mp.mpf(1)
    for i in range(dmax + 1):
        r = [0] * (dmax + 1)
        r[i] = 1
        r.append(int(mp.nint(C * pw)))
        rows.append(r)
        pw *= x
    M = flint.fmpz_mat(rows).lll()
    best = None
    for r in range(M.nrows()):
        co = [int(M[r, k]) for k in range(dmax + 1)]
        if all(c == 0 for c in co):
            continue
        poly = flint.fmpz_poly(co)
        val = sum(mp.mpf(c) * x ** k for k, c in enumerate(co))
        H = max(abs(c) for c in co)
        if abs(val) < mp.mpf(10) ** (-(prec_digits // 2)) * H:
            best = poly
            break
    if best is None:
        return None
    # irreducible factor vanishing at x
    fac = best.factor()[1]
    cand = []
    for f, e in fac:
        co = [int(c) for c in f.coeffs()]
        val = abs(sum(mp.mpf(c) * x ** k for k, c in enumerate(co)))
        cand.append((val, f))
    cand.sort(key=lambda a: a[0])
    return cand[0][1]


def galois_probe(f):
    co = [int(c) for c in f.coeffs()]
    deg = f.degree()
    lead = co[-1]
    disc = None
    patterns = {}
    for p in PRIMES:
        if lead % p == 0:
            continue
        g = flint.nmod_poly(co, p)
        if g.degree() != deg:
            continue
        fl = g.factor()[1]
        if any(e > 1 for _, e in fl):
            continue  # ramified / bad
        pat = tuple(sorted(h.degree() for h, _ in fl))
        patterns[pat] = patterns.get(pat, 0) + 1
    orders = set()
    for pat in patterns:
        l = 1
        for a in pat:
            l = l * a // math.gcd(l, a)
        orders.add(l)
    real_roots = sum(1 for r in f.complex_roots() if r[0].imag.contains(0) and abs(float(r[0].imag.mid())) < 1e-30)
    return dict(deg=deg, orders=sorted(orders), npat=len(patterns),
                equal_deg=all(len(set(p)) == 1 for p in patterns),
                real_roots=real_roots, poly=[int(c) for c in co] if deg <= 12 else None)


def work(o):
    try:
        T, res, step = refine([(float(a), float(b)) for a, b in o['coords']], o['cls'], o['n'])
        out = []
        for c, t in sorted(T.items()):
            if c == o['cls'][0]:
                continue
            f = minpoly(t, 16) or minpoly(t, 32)
            if f is None:
                out.append(dict(cls=c, value=mp.nstr(t, 20), deg=None))
                continue
            g = galois_probe(f)
            g.update(cls=c, value=mp.nstr(t, 20))
            out.append(g)
        return dict(cell=o['cell'], i=o['i'], label=o['label'], frame=o['frame'], miss=o['miss'],
                    res=res, step=step, ratios=out)
    except Exception as e:
        return dict(cell=o['cell'], i=o['i'], label=o['label'], err=repr(e))


if __name__ == '__main__':
    R = json.load(open(sys.argv[1]))
    if len(sys.argv) > 3 and sys.argv[3] == 'missing':
        R = [o for o in R if o['miss']]
    elif len(sys.argv) > 3 and sys.argv[3] == 'exact':
        R = [o for o in R if not o['miss']]
    with Pool(2) as P:
        res = []
        for k, r in enumerate(P.imap_unordered(work, R, chunksize=2)):
            res.append(r)
            if k % 25 == 0:
                print(k, len(R), flush=True)
                json.dump(res, open(sys.argv[2], 'w'))
    json.dump(res, open(sys.argv[2], 'w'))
    print('done', len(res))
