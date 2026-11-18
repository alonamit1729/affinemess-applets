"""Fourteen points in the plane with at most seven distances: all of them, up to similarity.

Every such set other than the regular 14-gon and the regular 15-gon minus a vertex contains
eight points with at most six distances (README.md, the diameter-graph lemma).  This program
takes every eight-point set with at most six distances and finds every way to add six points
without exceeding seven distances.

  * The 1565 eight-point patterns of certificates/L8.json are solved (fewdist.solve_lengths);
    the solutions are saved in solutions8.json.  1551 patterns have finitely many solutions;
    the other 14 have a one-parameter family of solutions, a conic in two squared lengths.
  * Each rigid realisation K is extended exactly: candidate points (all cases of the argument
    in `candidates`, including those in which the new point is irrational over the field of
    K), then a depth-first search for six candidates that fit together.
  * Each family K(s) is analysed with s symbolic, as in fifteen/fifteen.py: a fourteen-point
    set can contain K(s0) only for s0 among the roots of finitely many polynomials, and K(s0)
    is then extended exactly for each of them.

Result: 82 sets.  81 have squared distances proportional to 1, 3, 4, 7, 9, 12, 13, and they
are exactly the 81 sets with these squared distances (fewdist.sets_with_distances); the other
one is the regular 13-gon with its centre.  With the regular 14-gon and the regular 15-gon
minus a vertex: 84 sets.

    python3 fourteen.py [--jobs K] [--resolve]      (about an hour on two cores)
"""
import json
import multiprocessing as mp
import os
import sys
import time
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from math import isqrt

import mpmath
import sympy as sp
from mpmath import iv

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "fifteen"))
import fewdist as fd                                            # noqa: E402
import fifteen as ff                                            # noqa: E402
from fifteen import QS, _numerator, at_value                    # noqa: E402

NEED = 6                                                         # 8 + 6 = 14 points
PALETTE = [1, 3, 4, 7, 9, 12, 13]
iv.dps = 50


# ---- small geometry -------------------------------------------------------------------------

def rat_sqrt(c):
    """sqrt(c) if the rational number c is the square of a rational number, else None."""
    c = Fraction(c)
    if c < 0:
        return None
    a, b = isqrt(c.numerator), isqrt(c.denominator)
    return Fraction(a, b) if a * a == c.numerator and b * b == c.denominator else None


def line_row(conf, u, v, rho, sig):
    """The points y with |y-u|^2 - |y-v|^2 = rho - sig, as a row (a, b, c): a*y0 + b*y1 = c."""
    e = (v[0] - u[0], v[1] - u[1])
    z = (0, 0)
    return (2 * (conf.A * e[0] + conf.C * e[1]), 2 * (conf.C * e[0] + conf.B * e[1]),
            rho - sig + conf.d2(v, z) - conf.d2(u, z))


def normal(conf, e):
    """A vector orthogonal to e for the inner product of conf."""
    return (-(conf.C * e[0] + conf.B * e[1]), conf.A * e[0] + conf.C * e[1])


def anchors(conf, k):
    """Non-collinear triples such that every k points that are not collinear contain one of
    them (all non-collinear triples when k < 3)."""
    P, n = conf.points, len(conf.points)
    noncol = [T for T in combinations(range(n), 3) if not conf.collinear(*(P[i] for i in T))]
    if k < 3:
        return noncol
    mask = lambda T: sum(1 << i for i in T)
    masks = [mask(T) for T in noncol]
    todo = {mask(S) for S in combinations(range(n), k)
            if any(not conf.collinear(*(P[i] for i in T)) for T in combinations(S, 3))}
    chosen = []
    while todo:
        best = max(range(len(noncol)),
                   key=lambda j: sum((masks[j] & W) == masks[j] for W in todo))
        todo -= {W for W in todo if (masks[best] & W) == masks[best]}
        chosen.append(noncol[best])
    return chosen


def lines_of(conf):
    P, n = conf.points, len(conf.points)
    return {frozenset(i for i in range(n) if i in T or conf.collinear(P[T[0]], P[T[1]], P[i]))
            for T in combinations(range(n), 2)}


# ---- candidate points -----------------------------------------------------------------------

def candidates(conf, b):
    """Every point y that can be added to conf with at most b new squared distances.

    Returns (good, irrational): `good` maps the points with coordinates in the field of conf to
    their sets of new distances; `irrational` lists triples (f, nv, c) standing for the two
    points f +- sqrt(c) nv, where sqrt(c) may lie outside the field (decided later).

    Let S be the set of points of conf at old distances from y and r = n - |S|.
     (i)   If a value occurs three times among the n distances, y is a circumcentre.
     Otherwise every value occurs at most twice, so r <= 2b and at least r - b of the new
     values occur exactly twice.
     (ii)  If S contains three non-collinear points, then |S| >= n - 2b and S contains an anchor
           triple; y is determined by its old distances to it (two linear equations).
     (iii) Otherwise S is collinear (|S| >= n - 2b >= 2 since n = 8, b <= 3), which needs n - 2b
           collinear points in conf.  For a new value taken twice, y is on the perpendicular
           bisector of the two points.
        (a) Two such bisectors that meet determine y.
        (b1) If there is one such bisector L (all of them coincide; distinct parallel ones admit
             no y), take u, v in S: y also satisfies |y-u|^2 - |y-v|^2 = rho - sig, a line
             perpendicular to uv.  If it is not parallel to L, it meets L in y.  If it is, L is
             perpendicular to uv through the point f of the line uv, and y = f +- h nv with
             h^2 = rho - |f-u|^2, nv normal to uv.
        (b2) If no new value is taken twice, r <= b and S has n - b collinear points; with u, v
             in S, the foot f of y on their line is fixed by rho - sig, and y = f +- h nv as
             in (b1).
    """
    P = conf.points
    n = len(P)
    old = conf.distances()
    vals = sorted(old, key=repr)
    zero = vals[0] - vals[0]
    good, seen, irrational = {}, set(P), set()

    def consider(y):
        if y is None or y in seen:
            return
        seen.add(y)
        new = set()
        for p in P:
            d = conf.d2(y, p)
            if d not in old:
                new.add(d)
                if len(new) > b:
                    return
        good[y] = frozenset(new)

    def pair(f, nv, c):
        if c == 0:
            consider(f)
        elif isinstance(c, fd.Alg):
            irrational.add((f, nv, c))
        elif c > 0:
            s = rat_sqrt(c)
            if s is None:
                irrational.add((f, nv, c))
            else:
                consider((f[0] + s * nv[0], f[1] + s * nv[1]))
                consider((f[0] - s * nv[0], f[1] - s * nv[1]))

    noncol = [T for T in combinations(range(n), 3) if not conf.collinear(*(P[t] for t in T))]
    for T in noncol:                                                        # (i)
        consider(fd._locate(conf, T, zero, zero, zero))
    k = n - 2 * b
    for T in anchors(conf, k):                                              # (ii)
        for r, s, t in product(vals, repeat=3):
            y = fd._locate(conf, T, r, s, t)
            if conf.d2(y, P[T[0]]) == r:
                consider(y)
    lines = lines_of(conf)
    if max(len(L) for L in lines) < k:
        return good, irrational
    pairs = list(combinations(range(n), 2))
    bis = {pq: ff._bisector(conf, P[pq[0]], P[pq[1]]) for pq in pairs}
    for p1, p2 in combinations(pairs, 2):                                   # (iii a)
        if len(set(p1) | set(p2)) == 4:
            consider(ff._meet(bis[p1], bis[p2]))
    for ab in pairs:                                                        # (iii b1)
        L = bis[ab]
        for uv in pairs:
            if set(uv) & set(ab):
                continue
            u, v = P[uv[0]], P[uv[1]]
            if ff._meet(L, line_row(conf, u, v, vals[0], vals[0])) is not None:
                for rho, sig in product(vals, repeat=2):
                    consider(ff._meet(L, line_row(conf, u, v, rho, sig)))
            else:
                e = (v[0] - u[0], v[1] - u[1])
                a0, b0, c0 = L
                t = (c0 - a0 * u[0] - b0 * u[1]) / (a0 * e[0] + b0 * e[1])
                f = (u[0] + t * e[0], u[1] + t * e[1])
                nv = normal(conf, e)
                for rho in vals:
                    pair(f, nv, (rho - conf.d2(f, u)) / conf.dot(nv, nv))
    for M in lines:                                                         # (iii b2)
        if len(M) < n - b:
            continue
        i, j = sorted(M)[:2]
        u, v = P[i], P[j]
        e = (v[0] - u[0], v[1] - u[1])
        nv = normal(conf, e)
        for rho, sig in product(vals, repeat=2):
            rhs = rho - sig + conf.d2(v, (0, 0)) - conf.d2(u, (0, 0))
            t = (rhs - 2 * conf.dot(u, e)) / (2 * conf.dot(e, e))
            f = (u[0] + t * e[0], u[1] + t * e[1])
            pair(f, nv, (rho - conf.d2(f, u)) / conf.dot(nv, nv))
    return good, irrational


# ---- the points f +- sqrt(c) nv -------------------------------------------------------------

X_ = sp.Symbol("x")


def sqrt_in_field(c, F):
    """r in F with r^2 = c, or None if c is not a square in F.  For F = Q exactly; otherwise
    sympy factors x^2 - c over Q(alpha), and a square root found by PSLQ is checked exactly."""
    if not isinstance(c, fd.Alg):
        return rat_sqrt(c)
    a = sp.CRootOf(F.sympy_poly, 0)
    expr = sum(sp.Rational(x.numerator, x.denominator) * a ** k for k, x in enumerate(c.c))
    _, facs = sp.factor_list(X_ ** 2 - expr, X_, extension=a)
    if all(sp.degree(g, X_) == 2 for g, _ in facs):
        return None
    i = next(j for j in range(F.real_roots) if F.sign(c, j) > 0)
    mpmath.mp.dps = 60
    while F.intervals[i][1] - F.intervals[i][0] > Fraction(1, 10 ** 70):
        F._narrow(i)
    al = mpmath.mpf(F.intervals[i][0].numerator) / F.intervals[i][0].denominator
    cv = sum(mpmath.mpf(x.numerator) / x.denominator * al ** k for k, x in enumerate(c.c))
    rel = mpmath.pslq([mpmath.sqrt(cv)] + [al ** k for k in range(F.deg)],
                      maxcoeff=10 ** 12, maxsteps=10 ** 6)
    assert rel and rel[0] != 0, "PSLQ failed"
    r = F.element([Fraction(-rel[k + 1], rel[0]) for k in range(F.deg)])
    assert r * r == c
    return r


def irrational_ok(conf, b, item):
    """y = f + s nv with s^2 = c and s not in the field: does y have at most b new distances?
    |y-p|^2 = (|f-p|^2 + c<nv,nv>) + s * 2<nv,f-p> is compared exactly as a pair over the field
    (the conjugate point f - s nv gives the same answer)."""
    f, nv, c = item
    old = conf.distances()
    new = set()
    for p in conf.points:
        w = (f[0] - p[0], f[1] - p[1])
        a, bb = conf.dot(w, w) + c * conf.dot(nv, nv), 2 * conf.dot(nv, w)
        if not (bb == 0 and a in old):
            new.add((a, bb))
    return len(new) <= b


def to_iv(q):
    q = Fraction(q)
    return iv.mpf(q.numerator) / q.denominator


def enclose(a, F, i):
    """A certified interval for a at the i-th real root of the field F."""
    if not isinstance(a, fd.Alg):
        return to_iv(a)
    while F.intervals[i][1] - F.intervals[i][0] > Fraction(1, 10 ** 45):
        F._narrow(i)
    L, H = F.enclose(a, i)
    return iv.mpf([to_iv(L).a, to_iv(H).b])


def stab(intervals):
    """The least number of distinct reals that can lie in the given intervals, one in each."""
    count, last = 0, None
    for x in sorted(intervals, key=lambda x: float(x.b)):
        if last is None or x.a > last:
            count += 1
            last = x.b
    return count


def interval_hits(conf, b, field, root, good, irr):
    """Count the choices of six candidates, at least one of them irrational, that are not
    excluded by certified interval arithmetic (a choice is excluded only when the intervals of
    its new distances force more than b distinct values)."""
    E = (lambda a: enclose(a, field, root)) if field else to_iv
    A, B, C = E(conf.A), E(conf.B), E(conf.C)
    d2 = lambda p, q: (A * (p[0] - q[0]) ** 2 + 2 * C * (p[0] - q[0]) * (p[1] - q[1])
                       + B * (p[1] - q[1]) ** 2)
    P = [(E(x), E(y)) for x, y in conf.points]
    old = [E(v) for v in conf.distances()]
    pts = [((E(y[0]), E(y[1])), False) for y in good]
    for f, nv, c in irr:
        if field and field.sign(c, root) < 0:
            continue
        s = iv.sqrt(E(c))
        fx, nx = (E(f[0]), E(f[1])), (E(nv[0]), E(nv[1]))
        for sg in (1, -1):
            pts.append(((fx[0] + sg * s * nx[0], fx[1] + sg * s * nx[1]), True))
    if not any(t for _, t in pts):
        return 0
    newvals = lambda ds: [d for d in ds if all(d.b < o.a or d.a > o.b for o in old)]
    nv0 = [newvals([d2(y, p) for p in P]) for y, _ in pts]
    keep = [j for j in range(len(pts)) if stab(nv0[j]) <= b]
    hits = 0

    def search(start, chosen, vals):
        nonlocal hits
        if len(chosen) == NEED:
            hits += any(pts[j][1] for j in chosen)
            return
        for jj in range(start, len(keep)):
            j = keep[jj]
            v = vals + nv0[j] + newvals([d2(pts[j][0], pts[c][0]) for c in chosen])
            if stab(v) <= b:
                search(jj + 1, chosen + [j], v)

    search(0, [], [])
    return hits


# ---- extending one eight-point set ----------------------------------------------------------

def completions(conf, b, good, need=NEED):
    """All sets of `need` exact candidates that can be added together within the budget."""
    old = conf.distances()
    cand = list(good.items())
    found = []

    def search(start, chosen, new):
        if len(chosen) == need:
            found.append(tuple(chosen))
            return
        if len(cand) - start < need - len(chosen):
            return
        for j in range(start, len(cand)):
            y, new_y = cand[j]
            total = new | new_y
            if len(total) > b:
                continue
            for c in chosen:
                d = conf.d2(y, cand[c][0])
                if d not in old:
                    total = total | {d}
                    if len(total) > b:
                        break
            if len(total) <= b:
                search(j + 1, chosen + [j], total)

    search(0, [], frozenset())
    return [c[0] for c in cand], found


def extend(sol):
    """The shapes of all fourteen-point sets containing the eight-point solution `sol`, and the
    number of interval hits (choices involving irrational points not excluded; 0 throughout)."""
    roots = sol.embeddings
    if not roots:
        return set(), 0, 0
    K = fd.place(sol.pattern, sol.lengths)
    b = 7 - len(K.distances())
    F = sol.field
    good, irrational = candidates(K, b)
    irr = []
    for item in irrational:
        f, nv, c = item
        r = sqrt_in_field(c, F)
        if r is None:
            if irrational_ok(K, b, item):
                irr.append(item)
            continue
        for y in ((f[0] + r * nv[0], f[1] + r * nv[1]), (f[0] - r * nv[0], f[1] - r * nv[1])):
            if y not in good and y not in K.points:
                new = {K.d2(y, p) for p in K.points} - K.distances()
                if len(new) <= b:
                    good[y] = frozenset(new)
    cands, comps = completions(K, b, good)
    shapes = set()
    for root in roots:
        for c in comps:
            shapes.add(ff.shape_key(K.with_points(cands[j] for j in c), root, F))
    hits = sum(interval_hits(K, b, F, root, cands, irr) for root in roots) if irr else 0
    return shapes, hits, len(roots)


# ---- the eight-point patterns ---------------------------------------------------------------

def load_L8():
    return fd.load_patterns(os.path.join(HERE, "..", "certificates", "L8.json"))


def solve_all(patterns, path, resolve=False):
    """Solutions of the 1565 patterns; a flexible pattern is stored with its lex basis."""
    if os.path.exists(path) and not resolve:
        data = json.load(open(path))["by_pattern"]
    else:
        data = []
        for p in patterns:
            try:
                data.append({"solutions": [s.to_json() for s in fd.solve_lengths(p)]})
            except ValueError:
                xs, w, S = fd._saturated(p, "principal")
                G = sp.groebner(list(S.exprs), w, *xs, order="lex")
                data.append({"flexible": True,
                             "basis": [str(g) for g in G.exprs if not g.has(w)]})
        json.dump({"patterns": len(data), "by_pattern": data}, open(path, "w"),
                  separators=(",", ":"))
    fields, out = {}, []
    for p, d in zip(patterns, data):
        out.append(d["basis"] if d.get("flexible") else
                   [fd.LengthSolution.from_json(p, x, fields) for x in d["solutions"]])
    return out


# ---- the fourteen flexible patterns ---------------------------------------------------------

S_ = sp.Symbol("s")


def parametrise(pattern, basis):
    """The squared lengths of a flexible pattern as elements of Q(s), from its lex basis (linear
    equations and one conic g(p, q) = 0 in two of the lengths): lines of slope s through a
    rational point (p0, q0) of the conic.  Also returns the points of the conic that the
    parametrisation misses (those with p = p0), as lists of rational lengths."""
    xs, values = fd.length_symbols(pattern)
    X = {str(x): x for x in xs}
    G = [sp.sympify(g, locals=X) for g in basis]
    conic = [g for g in G if sp.Poly(g, *xs).total_degree() == 2]
    assert len(conic) == 1, basis
    g = conic[0]
    pv, qv = sorted(g.free_symbols, key=lambda v: xs.index(v))
    lin = [h for h in G if h is not g]
    assert all(sp.Poly(h, *xs).total_degree() == 1 for h in lin)
    sol_lin, = sp.solve(lin, [x for x in xs if x not in (pv, qv)], dict=True)
    base = None
    for den in range(1, 40):
        for num in range(-80, 81):
            q0 = sp.Rational(num, den)
            for r in sp.Poly(g.subs(qv, q0), pv).ground_roots():
                if r.is_rational:
                    base = (r, q0)
                    break
            if base:
                break
        if base:
            break
    assert base, f"no small rational point on {g}"
    p0, q0 = base
    t = sp.Symbol("t")
    h = sp.expand(g.subs({pv: p0 + t, qv: q0 + S_ * t}))
    tt = sp.cancel(-h.coeff(t, 1) / h.coeff(t, 2))

    def lengths_at(pp, qq):
        v = {pv: pp, qv: qq}
        for x, e in sol_lin.items():
            v[x] = sp.cancel(e.subs({pv: pp, qv: qq}))
        return v

    gen = lengths_at(sp.cancel(p0 + tt), sp.cancel(q0 + S_ * tt))
    for e in G:
        assert sp.simplify(e.subs(gen)) == 0
    conv = lambda e: QS.from_expr(sp.fraction(sp.cancel(e))[0]) / QS.from_expr(
        sp.fraction(sp.cancel(e))[1])
    lengths = [QS(1) if v == 1 else conv(gen[v]) for v in values]
    missed = []
    for r in sp.Poly(g.subs(pv, p0), qv).ground_roots():
        v = lengths_at(p0, r)
        missed.append([Fraction(1) if x == 1 else Fraction(int(sp.Rational(v[x]).p),
                                                           int(sp.Rational(v[x]).q))
                       for x in values])
    return lengths, missed, f"{g} = 0 through ({p0}, {q0})"


def _gcd_condition(newv, vals, b):
    """The gcd, over all (b+1)-subsets of the new values, of the product of the numerators of
    their differences and of their differences with the old values: it vanishes at every
    parameter value at which at most b new values remain."""
    g = None
    for group in combinations(newv, b + 1):
        q = None        # the numerator of a product divides the product of numerators
        for v1, v2 in list(combinations(group, 2)) + [(v1, x) for v1 in group for x in vals]:
            e = _numerator(v1 - v2)
            q = e if q is None else q * e
        g = q if g is None else g.gcd(q)
        if g.degree() == 0:
            break
    return g


def special_values(pattern, lengths, need=NEED):
    """Irreducible polynomials in s with real roots, whose roots include every s0 at which K(s0)
    can be part of a fourteen-point set (or at which the construction degenerates).  The
    candidate points of K(s) are built as in `candidates` with s symbolic.  Candidates that work
    for every s are allowed; then every `need` of them are checked together."""
    K = fd.place(pattern, lengths)
    P, n = K.points, len(K.points)
    old = K.distances()
    vals = sorted(old, key=str)
    b = 7 - len(old)
    zero = vals[0] - vals[0]
    polys = [_numerator(K.A * K.B - K.C * K.C)]
    noncol = [T for T in combinations(range(n), 3) if not K.collinear(*(P[t] for t in T))]
    for T in noncol:
        u, v, w = (P[t] for t in T)
        polys.append(_numerator((v[0] - u[0]) * (w[1] - u[1]) - (v[1] - u[1]) * (w[0] - u[0])))
    ys = [fd._locate(K, T, zero, zero, zero) for T in noncol]
    k = n - 2 * b
    for T in anchors(K, k):
        for r, t, u in product(vals, repeat=3):
            y = fd._locate(K, T, r, t, u)
            e = K.d2(y, P[T[0]]) - r
            if e == 0:
                ys.append(y)
            else:
                polys.append(_numerator(e))
    lines = lines_of(K)
    if max(len(L) for L in lines) >= k:
        pairs = list(combinations(range(n), 2))
        bis = {pq: ff._bisector(K, P[pq[0]], P[pq[1]]) for pq in pairs}
        det = lambda r1, r2: r1[0] * r2[1] - r1[1] * r2[0]
        for p1, p2 in combinations(pairs, 2):
            if len(set(p1) | set(p2)) == 4:
                d = det(bis[p1], bis[p2])
                if d != 0:
                    polys.append(_numerator(d))
                    ys.append(ff._meet(bis[p1], bis[p2]))
        for ab in pairs:
            for uv in pairs:
                if set(uv) & set(ab):
                    continue
                u, v = P[uv[0]], P[uv[1]]
                d = det(bis[ab], line_row(K, u, v, vals[0], vals[0]))
                assert d != 0, "case (iii b1) with s symbolic: not implemented"
                polys.append(_numerator(d))
                for rho, sig in product(vals, repeat=2):
                    ys.append(ff._meet(bis[ab], line_row(K, u, v, rho, sig)))
        assert all(len(M) < n - b for M in lines), "case (iii b2) with s symbolic"
    generic = []
    for y in set(ys):
        if y in P:
            continue
        newv = sorted({K.d2(y, q) for q in P} - old, key=str)
        if len(newv) <= b:
            generic.append(y)
        else:
            polys.append(_gcd_condition(newv, vals, b))
    for Y in combinations(generic, need):                  # candidates valid for every s
        newv = set()
        for y in Y:
            newv |= {K.d2(y, q) for q in P}
        for y, z in combinations(Y, 2):
            newv.add(K.d2(y, z))
        newv = sorted(newv - old, key=str)
        assert len(newv) > b, "six candidates that fit together for every s"
        polys.append(_gcd_condition(newv, vals, b))
    factors = set()
    for q in polys:
        if q is None or q.degree() <= 0:
            continue
        for f, _ in q.factor_list()[1]:
            if f.degree() > 0 and f.count_roots() > 0:
                factors.add(f.monic())
    return factors, len(generic)


# ---- main -----------------------------------------------------------------------------------

_WORK = {}


def _extend_index(args):
    i, k = args
    shapes, hits, cores = extend(_WORK["solutions"][i][k])
    return list(shapes), hits, cores


def _family(i):
    pattern, basis = _WORK["patterns"][i], _WORK["solutions"][i]
    lengths, missed, note = parametrise(pattern, basis)
    factors, generic = special_values(pattern, lengths)
    sols = [at_value(pattern, lengths, f) for f in factors]
    sols += [fd.LengthSolution(pattern, L) for L in missed
             if all(x != 0 for x in L) and len(set(L)) == len(L)]
    shapes, hits, members = set(), 0, 0
    for sol in sols:
        if sol is None or not sol.embeddings:
            continue
        members += 1
        sh, h, _ = extend(sol)
        shapes |= sh
        hits += h
    return i, note, generic, len(factors), members, list(shapes), hits


def regular13_with_centre_key():
    """The canonical named pattern of the regular 13-gon with its centre (ranks: chords 1, 2,
    then the radius, then chords 3 to 6)."""
    rank = {1: 0, 2: 1, 3: 3, 4: 4, 5: 5, 6: 6}
    pat = []
    for i, j in combinations(range(14), 2):
        if j == 13:
            pat.append(2)
        else:
            k = min(j - i, 13 - (j - i))
            pat.append(rank[k])
    return fd.canonical(tuple(pat), named=True)


def main():
    args = sys.argv[1:]
    jobs = int(args[args.index("--jobs") + 1]) if "--jobs" in args else 1
    t0 = time.time()
    patterns = load_L8()
    solutions = solve_all(patterns, os.path.join(HERE, "solutions8.json"),
                          "--resolve" in args)
    flexible = [i for i, x in enumerate(solutions) if x and isinstance(x[0], str)]
    print(f"{len(patterns)} eight-point patterns: {len(patterns) - len(flexible)} with "
          f"finitely many solutions, {len(flexible)} flexible: {flexible}  "
          f"({time.time() - t0:.0f} s)",
          flush=True)
    _WORK["solutions"], _WORK["patterns"] = solutions, patterns
    shapes, cores, hits = set(), 0, 0
    tasks = [(i, k) for i, x in enumerate(solutions) if i not in flexible
             for k in range(len(x))]
    fork = "fork" in mp.get_all_start_methods()
    ctx = mp.get_context("fork") if fork else mp.get_context()
    with ctx.Pool(jobs) as pool:
        for keys, h, c in pool.imap_unordered(_extend_index, tasks, chunksize=4):
            shapes.update(keys)
            hits += h
            cores += c
        print(f"rigid eight-point sets: {cores} realisations extended; {len(shapes)} "
              f"fourteen-point shapes; interval hits {hits}  ({time.time() - t0:.0f} s)",
              flush=True)
        for i, note, generic, nf, members, keys, h in pool.imap_unordered(_family, flexible):
            shapes.update(keys)
            hits += h
            print(f"family {i} ({note}): {generic} candidates for every s, {nf} special "
                  f"values, {members} genuine; {len(shapes)} shapes so far  "
                  f"({time.time() - t0:.0f} s)", flush=True)
    assert hits == 0
    target = tuple(("Q", Fraction(v)) for v in PALETTE)
    lattice = {k for k in shapes if k[0] == target}
    other = shapes - lattice
    found = fd.sets_with_distances(PALETTE, 14)
    assert {k[1] for k in lattice} == {fd.canonical(key, named=True) for key in found}
    assert len(other) == 1 and next(iter(other))[1] == regular13_with_centre_key()
    result = []
    for key, conf in found.items():
        result.append({"A": str(conf.A), "B": str(conf.B), "C": str(conf.C),
                       "points": [[str(x), str(y)] for x, y in conf.points],
                       "census": {str(k): v for k, v in sorted(Counter(
                           conf.d2(p, q) for p, q in combinations(conf.points, 2)).items())}})
    json.dump({"sets": result,
               "plus": ["the regular 13-gon with its centre", "the regular 14-gon",
                        "the regular 15-gon minus a vertex"]},
              open(os.path.join(HERE, "result.json"), "w"), indent=1)
    print(f"{len(shapes)} fourteen-point sets contain eight points with at most six distances: "
          f"{len(lattice)} with squared distances proportional to {PALETTE} (exactly the "
          f"{len(found)} sets with these squared distances; result.json) and the regular "
          f"13-gon with its centre.  With the regular 14-gon and the regular 15-gon minus a "
          f"vertex, "
          f"{len(shapes) + 2} sets  ({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
