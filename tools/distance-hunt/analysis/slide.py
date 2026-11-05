"""Can a catalog set be built with straightedge-and-compass steps plus 'slides':
put a point of the set on an existing line or circle, build what depends on it,
then equalize ONE pair of the set's distances (one equation pins the one parameter)?

Construction model (as in the game): curves through known points are lines through
two known points, perpendicular bisectors of two known points, and circles about a
known point whose radius is a distance between known points. A set point becomes
known when it lies on two different such curves (or when a slide pins it).
Every known point carries a recipe, so a slider can be perturbed and its dependents
recomputed; an equality between two of the set's distances pins the slider when it
is not identically true along the slide (nonzero derivative).

Usage: python3 slide.py rigid.json cross.json out.json
"""
import json, sys, math, cmath, itertools
from multiprocessing import Pool

TOL = 1e-8


def curve_val(spec, P):
    k = spec[0]
    if k == 'L':
        a, b = P[spec[1]], P[spec[2]]; d = b - a
        return ('l', a, d / abs(d))
    if k == 'B':
        a, b = P[spec[1]], P[spec[2]]; d = b - a
        return ('l', (a + b) / 2, 1j * d / abs(d))
    c, p, q = P[spec[1]], P[spec[2]], P[spec[3]]
    return ('c', c, abs(p - q))


def meet(A, B):
    if A[0] == 'l' and B[0] == 'l':
        (_, p, u), (_, q, v) = A, B
        den = (u.conjugate() * v).imag
        if abs(den) < 1e-12: return []
        t = ((q - p).conjugate() * v).imag / den
        return [p + t * u]
    if A[0] == 'c' and B[0] == 'l': A, B = B, A
    if A[0] == 'l':
        _, p, u = A; _, c, r = B
        t0 = ((c - p) * u.conjugate()).real; f = p + t0 * u
        h2 = r * r - abs(c - f) ** 2
        if h2 < -1e-12: return []
        h = math.sqrt(max(h2, 0)); return [f - h * u, f + h * u]
    _, c1, r1 = A; _, c2, r2 = B
    d = abs(c2 - c1)
    if d < 1e-12: return []
    a = (r1 * r1 - r2 * r2 + d * d) / (2 * d); h2 = r1 * r1 - a * a
    if h2 < -1e-12: return []
    e = (c2 - c1) / d; b = c1 + a * e; h = math.sqrt(max(h2, 0))
    return [b + 1j * h * e, b - 1j * h * e]


def on_curve(spec, P, z):
    v = curve_val(spec, P)
    if v[0] == 'l':
        _, p, u = v; return abs(((z - p) * u.conjugate()).imag) < TOL
    _, c, r = v; return abs(abs(z - c) - r) < TOL


def curve_key(spec, P):
    v = curve_val(spec, P)
    if v[0] == 'l':
        _, p, u = v
        if u.real < -1e-12 or (abs(u.real) < 1e-12 and u.imag < 0): u = -u
        off = ((p) * u.conjugate()).imag
        return ('l', round(u.real, 6), round(u.imag, 6), round(off, 6))
    _, c, r = v
    return ('c', round(c.real, 6), round(c.imag, 6), round(r, 6))


class Build:
    def __init__(self, cfg, fixed):
        self.cfg = cfg
        self.rec = {i: ('fix',) for i in fixed}     # recipes for known set points
        self.order = list(fixed)

    def copy(self):
        b = Build.__new__(Build); b.cfg = self.cfg; b.rec = dict(self.rec); b.order = list(self.order); return b

    def curves_through(self, z):
        """distinct curves (by geometry) through z built from known points"""
        K = self.order; P = self.cfg; out = {}
        for a in K:
            ra = abs(z - P[a])
            if ra < TOL: continue
            for p, q in itertools.combinations(K, 2):
                if abs(abs(P[p] - P[q]) - ra) < TOL:
                    s = ('C', a, p, q); out.setdefault(curve_key(s, P), s); break
        for a, b in itertools.combinations(K, 2):
            for s in (('L', a, b), ('B', a, b)):
                if on_curve(s, P, z): out.setdefault(curve_key(s, P), s)
        return list(out.values())

    def close(self):
        n = len(self.cfg); grew = True
        while grew:
            grew = False
            for i in range(n):
                if i in self.rec: continue
                cs = self.curves_through(self.cfg[i])
                if len(cs) >= 2:
                    # most transversal pair
                    best = None
                    for c1, c2 in itertools.combinations(cs, 2):
                        m = meet(curve_val(c1, self.cfg), curve_val(c2, self.cfg))
                        if len(m) == 2: sep = abs(m[0] - m[1])
                        elif len(m) == 1: sep = 1.0
                        else: continue
                        if best is None or sep > best[0]: best = (sep, c1, c2)
                    if best is None: continue
                    self.rec[i] = ('x', best[1], best[2]); self.order.append(i); grew = True
        return len(self.rec) == n

    def evaluate(self, slider, t):
        """positions of known points when the slider moves by t along its curve"""
        P = dict()
        for i in self.order:
            r = self.rec[i]
            if r[0] == 'fix': P[i] = self.cfg[i]
            elif r[0] == 'slide':
                v = curve_val(r[1], P)
                if v[0] == 'l': P[i] = self.cfg[i] + t * v[2]
                else:
                    _, c, rad = v; ang = cmath.phase(self.cfg[i] - c) + t / rad
                    P[i] = c + rad * cmath.exp(1j * ang)
            else:
                m = meet(curve_val(r[1], P), curve_val(r[2], P))
                if not m: return None
                P[i] = min(m, key=lambda z: abs(z - self.cfg[i]))
        return P


def pins(b, s):
    """equalities among the set's distances that pin slider s (nonzero derivative)"""
    h = 1e-6
    Pp, Pm = b.evaluate(s, h), b.evaluate(s, -h)
    if Pp is None or Pm is None: return []
    K = b.order; cfg = b.cfg
    pairs = list(itertools.combinations(K, 2))
    out = []
    for (p1, p2) in itertools.combinations(pairs, 2):
        if abs(abs(cfg[p1[0]] - cfg[p1[1]]) - abs(cfg[p2[0]] - cfg[p2[1]])) > TOL: continue
        dp = abs(Pp[p1[0]] - Pp[p1[1]]) - abs(Pp[p2[0]] - Pp[p2[1]])
        dm = abs(Pm[p1[0]] - Pm[p1[1]]) - abs(Pm[p2[0]] - Pm[p2[1]])
        if abs(dp - dm) / (2 * h) > 1e-4: out.append((p1, p2))
    return out


def search(b, depth, memo, steps):
    if b.close(): return steps
    if depth == 0: return None
    key = frozenset(b.rec)
    if key in memo: return None
    memo.add(key)
    n = len(b.cfg)
    for i in range(n):
        if i in b.rec: continue
        for spec in b.curves_through(b.cfg[i]):
            nb = b.copy(); nb.rec[i] = ('slide', spec); nb.order.append(i)
            nb.close()
            pp = pins(nb, i)
            if not pp: continue
            # the slider is pinned: from now on it is a known point
            r = search(nb, depth - 1, memo, steps + [(i, spec, pp[0], len(nb.rec))])
            if r is not None: return r
    return None


def starts(cfg):
    import chain
    n = len(cfg); out = []
    for kind in ('tri', 'sq'):
        for i, j in itertools.permutations(range(n), 2):
            h = chain.lattice_hits(cfg[i], cfg[j], cfg, kind)
            if len(h) >= 3: out.append(('lattice-' + kind, [k for k in range(n) if any(abs(cfg[k] - q) < 1e-7 for q in h)]))
    for m in range(3, 19):
        for i, j in itertools.combinations(range(n), 2):
            for k in range(1, m // 2 + 1):
                for sg in (1, -1):
                    c, Vv = chain.polygon(cfg[i], cfg[j], m, k, sg)
                    base = Vv + [c]
                    on = [q for q in range(n) if any(abs(cfg[q] - z) < 1e-7 for z in base)]
                    if len(on) >= 3: out.append((f'{m}-gon', on))
    for i, j in itertools.combinations(range(n), 2): out.append(('pair', [i, j]))
    seen = set(); res = []
    for name, f in out:
        k = (name, tuple(sorted(f)))
        if k not in seen: seen.add(k); res.append((name, sorted(f)))
    res.sort(key=lambda a: -len(a[1]))
    return res


def work(o):
    cfg = [complex(float(x), float(y)) for x, y in o['coords']]
    best = None
    for name, f in starts(cfg):
        for depth in (1, 2, 3):
            r = search(Build(cfg, f), depth, set(), [])
            if r is not None:
                if best is None or len(r) < best[1]: best = (name, len(r), r)
                break
        if best and best[1] == 1: break
    return dict(cell=o['cell'], i=o['i'], label=o['label'],
                ok=best is not None, start=best and best[0], slides=best and best[1],
                plan=best and [[s[0], list(s[1]), [list(s[2][0]), list(s[2][1])]] for s in best[2]])


if __name__ == '__main__':
    R = json.load(open(sys.argv[1])); X = json.load(open(sys.argv[2]))
    V = {(c, i): (v, k) for c, i, v, k, l in X}
    sel = [dict(o, verdict=V[(o['cell'], o['i'])][0]) for o in R if V[(o['cell'], o['i'])][1] == 'Z']
    if len(sys.argv) > 4: sel = sel[:int(sys.argv[4])]
    print(len(sel), 'sets', flush=True)
    with Pool(2) as Pl:
        res = []
        for k, r in enumerate(Pl.imap_unordered(work, sel, chunksize=1)):
            res.append(r)
            if k % 20 == 0: print(k, flush=True); json.dump(res, open(sys.argv[3], 'w'))
    json.dump(res, open(sys.argv[3], 'w'))
    print('done', sum(r['ok'] for r in res), 'of', len(res))
